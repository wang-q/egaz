#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Run;
use AlignDB::Util qw(:all);

use Time::Duration;
use File::Find::Rule;
use File::Spec;
use File::Basename;
use File::Remove qw(remove);
use Path::Class;
use String::Compare;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# executable file location
my $kent_bin  = "~/bin/x86_64";
my $phast_bin = "~/share/phast/bin";

# dirs
my ( $dir_target, $dir_query );
my $dir_lav = ".";

# .synNet.maf or .net.maf
my $syn;

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'            => \$help,
    'man'               => \$man,
    'bin|kent_bin=s'    => \$kent_bin,
    'phast|phast_bin=s' => \$phast_bin,
    'dt|dir_target=s'   => \$dir_target,
    'dq|dir_query=s'    => \$dir_query,
    'dl|dir_lav=s'      => \$dir_lav,
    'syn'               => \$syn,
    'p|parallel=i'      => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

my $dir_synnet = "$dir_lav/synNet";
my $dir_chain  = "$dir_lav/chain";
for ( $dir_synnet, $dir_chain ) {
    mkdir $_, 0777 unless -e $_;
}

my $t_prefix = basename($dir_target);
my $q_prefix = basename($dir_query);

#----------------------------------------------------------#
# maf section
#----------------------------------------------------------#
if ( !$syn ) {
    my $dir_mafnet = "$dir_lav/mafNet";
    mkdir $dir_mafnet, 0777 unless -e $dir_mafnet;

    #----------------------------#
    # axtToMaf section
    #----------------------------#

    my @files = File::Find::Rule->file->name('*.axt')->in($dir_lav);
    @files = File::Find::Rule->file->name('*.axt.gz')->in($dir_lav)
        if scalar @files == 0;
    printf "\n----%4s .axt files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file = $job;
        my $output = basename( $file, ".axt", ".axt.gz" );
        $output = File::Spec->catfile( $dir_mafnet, "$output.maf" );

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
            = "$kent_bin/axtToMaf"
            . " -tPrefix=$t_prefix."
            . " -qPrefix=$q_prefix."
            . " $file"
            . " $dir_target/chr.sizes"
            . " $dir_query/chr.sizes"
            . " $output";
        exec_cmd($cmd);
        print ".maf file generated.\n\n";

        return;
    };

    my @jobs = sort @files;
    my $run  = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
}
else {
    my $dir_mafsynnet = "$dir_lav/mafSynNet";
    mkdir $dir_mafsynnet, 0777 unless -e $dir_mafsynnet;

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
        my $cmd
            = "$kent_bin/netFilter" . " -syn"
            . " $files[0]"
            . " | $kent_bin/netSplit stdin"
            . " $dir_synnet";
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
        my $cmd = "$kent_bin/chainSplit" . " $dir_chain" . " $files[0]";
        exec_cmd($cmd);
    }

    my @files = File::Find::Rule->file->name('*.net')->in($dir_synnet);
    printf "\n----%4s .net files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file   = $job;
        my $base   = basename( $file, ".net" );
        my $output = File::Spec->catfile( $dir_mafsynnet, "$base.synNet.maf" );
        my $chain_file = File::Spec->catfile( $dir_chain, "$base.chain" );

        print "Run netToAxt axtSort axtToMaf...\n";
        my $cmd
            = "$kent_bin/netToAxt" 
            . " $file"
            . " $chain_file"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " stdout"
            . " | $kent_bin/axtSort stdin stdout"
            . " | $kent_bin/axtToMaf"
            . " -tPrefix=$t_prefix."
            . " -qPrefix=$q_prefix."
            . " stdin"
            . " $dir_target/chr.sizes"
            . " $dir_query/chr.sizes"
            . " $output";
        exec_cmd($cmd);
        print ".maf file generated.\n\n";

        return;
    };

    my @jobs = sort @files;
    my $run  = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
}

#----------------------------------------------------------#
# clean section
#----------------------------------------------------------#
{
    chdir $dir_lav;
    remove( \1, $dir_synnet, $dir_chain );

    my @files = File::Find::Rule->file->name('*.maf')->in("$dir_lav");
    printf "\n----%4s .maf files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;

        my $file = $job;
        my $cmd  = "gzip $file";
        exec_cmd($cmd);

        return;
    };

    my @jobs = sort @files;
    my $run  = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;
}

print "\n";
print "Runtime ", duration( time - $start_time ), ".\n";
print "=" x 30, "\n";

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

__END__

=head1 NAME

    amp.pl - axt-maf-phast pipeline

=head1 SYNOPSIS

    amp.pl -dt <one target dir or file> -dq <one query dir or file> [options]
      Options:
        -?, --help              brief help message
        --man                   full documentation

      Run in parallel mode
        -p, --paralle           number of child processes

      Fasta dirs  
        -dt, --dir_target       dir of target fasta files
        -dq, --dir_query        dir of query fasta files

      Output .lav and .axt
        -dl, --dir_lav          where .lav and .axt files storess

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION


=cut
