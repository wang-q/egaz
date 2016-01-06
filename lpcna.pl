#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;

use File::Find::Rule;
use File::Remove qw(remove);
use Path::Tiny;
use Time::Duration;
use IPC::Cmd qw(can_run);

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

lpcna.pl - lav-psl-chain-net-axt pipeline

=head1 SYNOPSIS

    perl lpcna.pl -dt <target> -dq <query> [options]
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
        --linearGap         STR     axtChain linearGap, loose or medium, default is [loose]
                                    Human18vsChimp2 use loose and 1000
                                    Human19vsChimp3 use medium and 5000
                                    loose is chicken/human linear gap costs.
                                    medium is mouse/human linear gap costs.
                                    Or specify a piecewise linearGap tab delimited file.
        --minScore          INT     Minimum score for axtChain

=cut

GetOptions(
    'help|?'          => sub { HelpMessage(0) },
    'dir_target|dt=s' => \my $dir_target,
    'dir_query|dq=s'  => \my $dir_query,
    'dir_lav|dl=s' => \( my $dir_lav   = '.' ),
    'parallel|p=i' => \( my $parallel  = 1 ),
    'linearGap=i'  => \( my $linearGap = "loose" ),
    'minScore=i'   => \( my $minScore  = "1000" ),
) or HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my @executables = qw{
    faops faToTwoBit lavToPsl axtChain chainAntiRepeat chainMergeSort chainPreNet
    chainNet netSyntenic netChainSubset chainStitchId netSplit netToAxt axtSort
};

for (@executables) {
    if ( !can_run($_) ) {
        die "Can't find $_\n";
    }
}

my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

# make lav dir
my $dir_net    = "$dir_lav/net";
my $dir_axtnet = "$dir_lav/axtNet";
for ( $dir_lav, $dir_net, $dir_axtnet ) {
    $_ = path($_)->absolute->stringify;
    path($_)->mkpath;
}

#----------------------------------------------------------#
# fa section
#----------------------------------------------------------#
{

    # faops size
    my $cmd = "faops size" . " $dir_target/*.fa" . " > $dir_target/chr.sizes";
    exec_cmd($cmd) if !-e "$dir_target/chr.sizes";

    $cmd = "faops size" . " $dir_query/*.fa" . " > $dir_query/chr.sizes";
    exec_cmd($cmd) if !-e "$dir_query/chr.sizes";

    # use combined .2bit file instead of dir of nibs

    # faToTwoBit - Convert DNA from fasta to 2bit format
    # usage:
    #    faToTwoBit in.fa [in2.fa in3.fa ...] out.2bit
    $cmd = "faToTwoBit" . " $dir_target/*.fa" . " $dir_target/chr.2bit";
    exec_cmd($cmd) if !-e "$dir_target/chr.2bit";

    $cmd = "faToTwoBit" . " $dir_query/*.fa" . " $dir_query/chr.2bit";
    exec_cmd($cmd) if !-e "$dir_query/chr.2bit";

    print "\n";
}

#----------------------------------------------------------#
# lavToPsl section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.lav')->in($dir_lav);
    printf "\n----%4s .lav files to be converted ----\n", scalar @files;

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file   = $chunk_ref->[0];
            my $output = $file;
            $output =~ s/lav$/psl/;

            # lavToPsl - Convert lav to psl format
            # usage:
            #   lavToPsl in.lav out.psl
            print "Run lavToPsl...\n";
            my $cmd = "lavToPsl" . " $file" . " $output";
            exec_cmd($cmd);
        }
    );

    print "\n";
}

#----------------------------------------------------------#
# axtChain section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.psl')->in($dir_lav);
    printf "\n----%4s .psl files to be converted ----\n", scalar @files;

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file   = $chunk_ref->[0];
            my $output = $file;
            $output =~ s/psl$/chain/;

            # axtChain - Chain together axt alignments.
            # usage:
            #   axtChain -linearGap=loose in.axt tNibDir qNibDir out.chain
            # Where tNibDir/qNibDir are either directories full of nib files, or the
            # name of a .2bit file
            #
            # chainAntiRepeat - Get rid of chains that are primarily the results of
            # repeats and degenerate DNA
            # usage:
            #    chainAntiRepeat tNibDir qNibDir inChain outChain
            # options:
            #    -minScore=N - minimum score (after repeat stuff) to pass
            #    -noCheckScore=N - score that will pass without checks (speed tweak)
            print "Run axtChain...\n";
            my $cmd
                = "axtChain -minScore=$minScore -linearGap=$linearGap -psl"
                . " $file"
                . " $dir_target/chr.2bit"
                . " $dir_query/chr.2bit"
                . " stdout"
                . " | chainAntiRepeat"
                . " $dir_target/chr.2bit"
                . " $dir_query/chr.2bit"
                . " stdin"
                . " $output";
            exec_cmd($cmd);
            print ".chain file generated.\n\n";
        }
    );

    print "\n";

    # This step would open all .chain files and reach system's maxfile limit.
    # So merge 100 files a time.
    #
    # chainMergeSort - Combine sorted files into larger sorted file
    # usage:
    #    chainMergeSort file(s)
    # Output goes to standard output
    # options:
    #    -saveId - keep the existing chain ids.
    #    -inputList=somefile - somefile contains list of input chain files.
    #    -tempDir=somedir/ - somedir has space for temporary sorting data, default ./
    my @files_chain = File::Find::Rule->file->name('*.chain')->in($dir_lav);
    my $i           = 1;
    while ( scalar @files_chain ) {
        my @batching = splice @files_chain, 0, 100;

        path( $dir_lav, "chainlist.tmp" )->spew( [ map {"$_\n"} @batching ] );

        my $cmd
            = "chainMergeSort"
            . " -inputList="
            . path( $dir_lav, "chainlist.tmp" )->stringify
            . " > $dir_lav/all.$i.chain.tmp";
        exec_cmd($cmd);
        path( $dir_lav, "chainlist.tmp" )->remove;

        $i++;
    }
    my $cmd = "chainMergeSort" . " $dir_lav/all.*.chain.tmp" . " > $dir_lav/all.chain";
    exec_cmd($cmd);

    {
        # chainPreNet - Remove chains that don't have a chance of being netted
        # usage:
        #   chainPreNet in.chain target.sizes query.sizes out.chain
        my $cmd
            = "chainPreNet"
            . " $dir_lav/all.chain"
            . " $dir_target/chr.sizes"
            . " $dir_query/chr.sizes"
            . " $dir_lav/all.pre.chain";
        exec_cmd($cmd);
    }

    print "\n";
}

#----------------------------------------------------------#
# chain-net section
#----------------------------------------------------------#
{

    # chainNet - Make alignment nets out of chains
    # usage:
    #   chainNet in.chain target.sizes query.sizes target.net query.net
    #
    # netSyntenic - Add synteny info to net.
    # usage:
    #   netSyntenic in.net out.net
    my $cmd
        = "chainNet -minSpace=1"
        . " $dir_lav/all.pre.chain"
        . " $dir_target/chr.sizes"
        . " $dir_query/chr.sizes"
        . " stdout"    # $dir_lav/target.chainnet
        . " /dev/null" . " | netSyntenic" . " stdin" . " $dir_lav/noClass.net";
    exec_cmd($cmd);

    # netChainSubset - Create chain file with subset of chains that appear in
    # the net
    # usage:
    #    netChainSubset in.net in.chain out.chain
    # options:
    #    -gapOut=gap.tab - Output gap sizes to file
    #    -type=XXX - Restrict output to particular type in net file
    #    -splitOnInsert - Split chain when get an insertion of another chain
    #    -wholeChains - Write entire chain references by net, don't split
    #     when a high-level net is encoundered.  This is useful when nets
    #     have been filtered.
    #    -skipMissing - skip chains that are not found instead of generating
    #     an error.  Useful if chains have been filtered.
    #
    # chainStitchId - Join chain fragments with the same chain ID into a single
    #    chain per ID.  Chain fragments must be from same original chain but
    #    must not overlap.  Chain fragment scores are summed.
    # usage:
    #    chainStitchId in.chain out.chain
    $cmd
        = "netChainSubset -verbose=0 $dir_lav/noClass.net"
        . " $dir_lav/all.chain"
        . " stdout"
        . " | chainStitchId"
        . " stdin"
        . " $dir_lav/over.chain";
    exec_cmd($cmd);

    # netSplit - Split a genome net file into chromosome net files
    # usage:
    #   netSplit in.net outDir
    $cmd = "netSplit" . " $dir_lav/noClass.net" . " $dir_net";
    exec_cmd($cmd);

    print "\n";
}

#----------------------------------------------------------#
# netToAxt section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.net')->in($dir_net);
    printf "\n----%4s .net files to be converted ----\n", scalar @files;

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file   = $chunk_ref->[0];
            my $output = path($file)->basename;
            $output .= ".axt";

            # netToAxt - Convert net (and chain) to axt.
            # usage:
            #   netToAxt in.net in.chain target.2bit query.2bit out.axt
            # note:
            # directories full of .nib files (an older format)
            # may also be used in place of target.2bit and query.2bit.
            #
            # axtSort - Sort axt files
            # usage:
            #   axtSort in.axt out.axt
            print "Run netToAxt...\n";
            my $cmd
                = "netToAxt"
                . " $file"
                . " $dir_lav/all.pre.chain"
                . " $dir_target/chr.2bit"
                . " $dir_query/chr.2bit"
                . " stdout"
                . " | axtSort stdin"
                . " $dir_axtnet/$output";
            exec_cmd($cmd);
            print ".axt file generated.\n\n";
        }
    );

    print "\n";
}

#----------------------------------------------------------#
# clean section
#----------------------------------------------------------#
{
    chdir $dir_lav;
    my $cmd;

    $cmd = "tar -czvf lav.tar.gz *.lav";    # bsdtar (mac) doesn't support  --remove-files
    if ( !-e "$dir_lav/lav.tar.gz" ) {
        exec_cmd($cmd);
        remove("*.lav");
    }

    remove( \1, "$dir_lav/net" );
    remove("[*.psl");
    remove("[*.chain");
    remove("*.tmp");

    $cmd = "gzip *.chain";
    exec_cmd($cmd);

    my @files = File::Find::Rule->file->name('*.axt')->in($dir_lav);
    printf "\n----%4s .axt files to be gzipped ----\n", scalar @files;

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file = $chunk_ref->[0];

            my $cmd = "gzip $file";
            exec_cmd($cmd);
        }
    );
}

print "\n";
print "Runtime ", duration( time - $start_time ), ".\n";
print "=" x 30, "\n";

exit;

__END__
