#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Run;
use AlignDB::Util qw(:all);

use Time::Duration;
use File::Find::Rule;
use File::Basename;
use File::Remove qw(remove);
use Path::Class;
use String::Compare;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# executable file location
my $kent_bin = "~/bin/x86_64";

# dirs
my ( $dir_target, $dir_query );
my $dir_lav = ".";

# Human18vsChimp2 use loose and 1000
# Human19vsChimp3 use medium and 5000
# axtChain linearGap
# loose is chicken/human linear gap costs.
# medium is mouse/human linear gap costs.
# Or specify a piecewise linearGap tab delimited file.
my $linearGap = "loose";

# axtChain minScore
# Minimum score for chain
my $minScore = "1000";

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'bin|kent_bin=s'  => \$kent_bin,
    'dt|dir_target=s' => \$dir_target,
    'dq|dir_query=s'  => \$dir_query,
    'dl|dir_lav=s'    => \$dir_lav,
    'l|linearGap=i'   => \$linearGap,
    'm|minScore=i'    => \$minScore,
    'p|parallel=i'    => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

# make dirs
unless ( -e $dir_lav ) {
    mkdir $dir_lav, 0777
        or die "Cannot create \"$dir_lav\" directory: $!";
}

my $dir_net    = "$dir_lav/net";
my $dir_axtnet = "$dir_lav/axtNet";
for ( $dir_net, $dir_axtnet ) {
    unless ( -e $_ ) {
        mkdir $_, 0777
            or die "Cannot create \"$_\" directory: $!";
    }
}

#----------------------------------------------------------#
# fa section
#----------------------------------------------------------#
{

    # faSize - print total base count in fa files.
    # usage:
    #   faSize file(s).fa
    my $cmd
        = "$kent_bin/faSize -detailed"
        . " $dir_target/*.fa"
        . " > $dir_target/chr.sizes";
    exec_cmd($cmd) if !-e "$dir_target/chr.sizes";

    $cmd
        = "$kent_bin/faSize -detailed"
        . " $dir_query/*.fa"
        . " > $dir_query/chr.sizes";
    exec_cmd($cmd) if !-e "$dir_query/chr.sizes";

    # use combined .2bit file instead of dir of nibs

    # faToTwoBit - Convert DNA from fasta to 2bit format
    # usage:
    #    faToTwoBit in.fa [in2.fa in3.fa ...] out.2bit
    $cmd
        = "$kent_bin/faToTwoBit"
        . " $dir_target/*.fa"
        . " $dir_target/chr.2bit";
    exec_cmd($cmd) if !-e "$dir_target/chr.2bit";

    $cmd = "$kent_bin/faToTwoBit" . " $dir_query/*.fa" . " $dir_query/chr.2bit";
    exec_cmd($cmd) if !-e "$dir_query/chr.2bit";

    print "\n";
}

#----------------------------------------------------------#
# lavToPsl section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.lav')->in($dir_lav);
    printf "\n----%4s .lav files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file   = $job;
        my $output = $file;
        $output =~ s/lav$/psl/;

        # lavToPsl - Convert blastz lav to psl format
        # usage:
        #   lavToPsl in.lav out.psl
        print "Run lavToPsl...\n";
        my $cmd = "$kent_bin/lavToPsl" . " $file" . " $output";
        exec_cmd($cmd);
        print ".psl file generated.\n\n";

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
# axtChain section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.psl')->in($dir_lav);
    printf "\n----%4s .psl files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file   = $job;
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
            = "$kent_bin/axtChain -minScore=$minScore -linearGap=$linearGap -psl"
            . " $file"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " stdout"
            . " | $kent_bin/chainAntiRepeat"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " stdin"
            . " $output";
        exec_cmd($cmd);
        print ".chain file generated.\n\n";

        return;
    };

    my @jobs = sort @files;

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
}

#----------------------------------------------------------#
# chain-net section
#----------------------------------------------------------#
{

    # chainMergeSort - Combine sorted files into larger sorted file
    # usage:
    #   chainMergeSort file(s)
    my $cmd
        = "$kent_bin/chainMergeSort"
        . " $dir_lav/*.chain"
        . " > $dir_lav/all.chain";
    exec_cmd($cmd);

    # chainPreNet - Remove chains that don't have a chance of being netted
    # usage:
    #   chainPreNet in.chain target.sizes query.sizes out.chain
    $cmd
        = "$kent_bin/chainPreNet"
        . " $dir_lav/all.chain"
        . " $dir_target/chr.sizes"
        . " $dir_query/chr.sizes"
        . " $dir_lav/all.pre.chain";
    exec_cmd($cmd);

    # chainNet - Make alignment nets out of chains
    # usage:
    #   chainNet in.chain target.sizes query.sizes target.net query.net
    #
    # netSyntenic - Add synteny info to net.
    # usage:
    #   netSyntenic in.net out.net
    $cmd
        = "$kent_bin/chainNet -minSpace=1"
        . " $dir_lav/all.pre.chain"
        . " $dir_target/chr.sizes"
        . " $dir_query/chr.sizes"
        . " stdout"    # $dir_lav/target.chainnet
        . " /dev/null"
        . " | $kent_bin/netSyntenic"
        . " stdin"
        . " $dir_lav/noClass.net";
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
        = "$kent_bin/netChainSubset -verbose=0 $dir_lav/noClass.net"
        . " $dir_lav/all.chain"
        . " stdout"
        . " | $kent_bin/chainStitchId"
        . " stdin"
        . " $dir_lav/over.chain";
    exec_cmd($cmd);

    # netSplit - Split a genome net file into chromosome net files
    # usage:
    #   netSplit in.net outDir
    $cmd = "$kent_bin/netSplit" . " $dir_lav/noClass.net" . " $dir_net";
    exec_cmd($cmd);

    print "\n";
}

#----------------------------------------------------------#
# netToAxt section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.net')->in($dir_net);
    printf "\n----%4s .net files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file   = $job;
        my $output = basename($file);
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
            = "$kent_bin/netToAxt" 
            . " $file"
            . " $dir_lav/all.pre.chain"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " stdout"
            . " | $kent_bin/axtSort stdin"
            . " $dir_axtnet/$output";
        exec_cmd($cmd);
        print ".axt file generated.\n\n";

        return;
    };

    my @jobs = sort @files;

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

}

#----------------------------------------------------------#
# clean section
#----------------------------------------------------------#
{
    chdir $dir_lav;
    my $cmd;

    $cmd = "tar -czvf lav.tar.gz [*.lav --remove-files";
    exec_cmd($cmd);

    remove( \1, "$dir_lav/net" );
    remove("[*.psl");
    remove("[*.chain");

    $cmd = "gzip *.chain";
    exec_cmd($cmd);

    my @files = File::Find::Rule->file->name('*.axt')->in("$dir_lav/axtNet");
    printf "\n----%4s .axt files to be converted ----\n", scalar @files;

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

    lpcna.pl - lav-psl-chain-net-axt pipeline

=head1 SYNOPSIS

    lpcna.pl -dt <one target dir or file> -dq <one query dir or file> [options]
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

-p 8
Runtime 3 minutes and 24 seconds.

-p 4
Runtime 3 minutes and 28 seconds.

-p 2
Runtime 4 minutes and 28 seconds.

-p 1
Runtime 8 minutes and 25 seconds.

=cut
