#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;

use File::Find::Rule;
use Path::Tiny;
use Time::Duration;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

amp.pl - axt-maf-phast pipeline

=head1 SYNOPSIS

    perl amp.pl -dt <target> -dq <query> [options]
      Options:
        --help          -?          brief help message

      Running mode
        --parallel      -p  INT     run in parallel mode, [1]

      Fasta dirs  
        --dir_target    -dt STR     dir of target fasta files
        --dir_query     -dq STR     dir of query fasta files

      Output .lav and .axt
        --dir_lav       -dl STR     where .lav and .axt files stores

      Chaining options
        --syn                       .synNet.maf or .net.maf

=cut

GetOptions(
    'help|?'          => sub { HelpMessage(0) },
    'dir_target|dt=s' => \my $dir_target,
    'dir_query|dq=s'  => \my $dir_query,
    'dir_lav|dl=s' => \( my $dir_lav  = '.' ),
    'parallel|p=i' => \( my $parallel = 1 ),
    'syn'          => \my $syn,
) or HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

my $dir_synnet = "$dir_lav/synNet";
my $dir_chain  = "$dir_lav/chain";
for ( $dir_synnet, $dir_chain ) {
    path($_)->mkpath;
}

my $t_prefix = path($dir_target)->basename;
my $q_prefix = path($dir_query)->basename;

#----------------------------------------------------------#
# maf section
#----------------------------------------------------------#
if ( !$syn ) {
    my $dir_mafnet = "$dir_lav/mafNet";
    path($dir_mafnet)->mkpath;

    #----------------------------#
    # axtToMaf section
    #----------------------------#

    my @files = File::Find::Rule->file->name('*.axt')->in($dir_lav);
    @files = File::Find::Rule->file->name('*.axt.gz')->in($dir_lav)
        if scalar @files == 0;
    printf "\n----%4s .axt files to be converted ----\n", scalar @files;

    my $worker = sub {
        my ( $self, $chunk_ref, $chunk_id ) = @_;

        my $file = $chunk_ref->[0];
        my $output = path($file)->basename( ".axt", ".axt.gz" );
        $output = path( $dir_mafnet, "$output.maf.gz" )->stringify;

        # axtToMaf - Convert from axt to maf format
        # usage:
        #    axtToMaf in.axt tSizes qSizes out.maf
        # Where tSizes and qSizes is a file that contains
        # the sizes of the target and query sequences.
        # Very often this with be a chrom.sizes file
        # Options:
        #     -qPrefix=XX. - add XX. to start of query sequence name in maf
        #     -tPrefix=YY. - add YY. to start of target sequence name in maf
        #     -tSplit Create a separate maf file for each target sequence.
        #             In this case output is a dir rather than a file
        #             In this case in.maf must be sorted by target.
        #     -score       - recalculate score
        #     -scoreZero   - recalculate score if zero
        print "Run axtToMaf...\n";
        my $cmd
            = "axtToMaf"
            . " -tPrefix=$t_prefix."
            . " -qPrefix=$q_prefix."
            . " $file"
            . " $dir_target/chr.sizes"
            . " $dir_query/chr.sizes"
            . " stdout"
            . " | gzip -c >"
            . " $output";
        exec_cmd($cmd);
        print ".maf file generated.\n\n";
    };

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach( [ sort @files ], $worker );

    print "\n";
}
else {
    my $dir_mafsynnet = "$dir_lav/mafSynNet";
    path($dir_mafsynnet)->mkpath;

    #----------------------------#
    # synNetMaf section
    #----------------------------#

    {
        my @files = File::Find::Rule->file->name('noClass.net')->in($dir_lav);
        @files = File::Find::Rule->file->name('noClass.net.gz')->in($dir_lav)
            if scalar @files == 0;
        if ( scalar @files == 0 ) {
            die "Can't find noClass.net\n";
        }

        # netFilter - Filter out parts of net.  What passes
        # filter goes to standard output.  Note a net is a
        # recursive data structure.  If a parent fails to pass
        # the filter, the children are not even considered.
        # usage:
        #    netFilter in.net(s)
        my $cmd = "netFilter" . " -syn" . " $files[0]" . " | netSplit stdin" . " $dir_synnet";
        exec_cmd($cmd);
    }

    {
        my @files = File::Find::Rule->file->name('all.chain')->in($dir_lav);
        @files = File::Find::Rule->file->name('all.chain.gz')->in($dir_lav)
            if scalar @files == 0;
        if ( scalar @files == 0 ) {
            die "Can't find all.chain\n";
        }

        # chainSplit - Split chains up by target or query sequence
        # usage:
        #    chainSplit outDir inChain(s)
        # options:
        #    -q  - Split on query (default is on target)
        #    -lump=N  Lump together so have only N split files.
        my $cmd = "chainSplit" . " $dir_chain" . " $files[0]";
        exec_cmd($cmd);
    }

    my @files = File::Find::Rule->file->name('*.net')->in($dir_synnet);
    printf "\n----%4s .net files to be converted ----\n", scalar @files;

    my $worker = sub {
        my ( $self, $chunk_ref, $chunk_id ) = @_;

        my $file       = $chunk_ref->[0];
        my $base       = path($file)->basename(".net");
        my $output     = path( $dir_mafsynnet, "$base.synNet.maf.gz" )->stringify;
        my $chain_file = path( $dir_chain, "$base.chain" )->stringify;

        print "Run netToAxt axtSort axtToMaf...\n";
        my $cmd
            = "netToAxt"
            . " $file"
            . " $chain_file"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " stdout"
            . " | axtSort stdin stdout"
            . " | axtToMaf"
            . " -tPrefix=$t_prefix."
            . " -qPrefix=$q_prefix."
            . " stdin"
            . " $dir_target/chr.sizes"
            . " $dir_query/chr.sizes"
            . " stdout"
            . " | gzip -c >"
            . " $output";
        exec_cmd($cmd);
        print ".maf file generated.\n\n";

        return;
    };

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach( [ sort @files ], $worker );

    print "\n";
}

#----------------------------------------------------------#
# clean section
#----------------------------------------------------------#
{
    path($dir_synnet)->remove_tree;
    path($dir_chain)->remove_tree;

}

print "\n";
print "Runtime ", duration( time - $start_time ), ".\n";
print "=" x 30, "\n";

exit;

__END__
