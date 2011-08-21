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

# axtChain linearGap
# loose is chicken/human linear gap costs.
# medium is mouse/human linear gap costs.
# Or specify a piecewise linearGap tab delimited file.
my $linearGap = "medium";

# axtChain minScore
# Minimum score for chain
my $minScore = "5000";

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'k|kent_bin=s'    => \$kent_bin,
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
# make lav dir
unless ( -e $dir_lav ) {
    mkdir $dir_lav, 0777
        or die "Cannot create \"$dir_lav\" directory: $!";
}

my $dir_per_chr = "$dir_lav/chr";

#----------------------------------------------------------#
# fa section
#----------------------------------------------------------#
{
    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    # faSize - print total base count in fa files.
    # usage:
    #   faSize file(s).fa
    my $cmd
        = "$kent_bin/faSize -detailed"
        . " $dir_target/*.fa"
        . " > $dir_target/chr.sizes";
    exec_cmd($cmd);

    $cmd
        = "$kent_bin/faSize -detailed"
        . " $dir_query/*.fa"
        . " > $dir_query/chr.sizes";
    exec_cmd($cmd);

    # faToTwoBit - Convert DNA from fasta to 2bit format
    # usage:
    #    faToTwoBit in.fa [in2.fa in3.fa ...] out.2bit
    $cmd
        = "$kent_bin/faToTwoBit"
        . " $dir_target/*.fa"
        . " $dir_target/chr.2bit";
    exec_cmd($cmd);
    
    $cmd
        = "$kent_bin/faToTwoBit"
        . " $dir_query/*.fa"
        . " $dir_query/chr.2bit";
    exec_cmd($cmd);

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

# use combined .2bit file instead of dir of nibs
#{
#    my @files
#        = File::Find::Rule->file->name('*.fa')->in( $dir_target, $dir_query );
#    printf "\n----%4s .fa files to be converted ----\n", scalar @files;
#
#    my $worker = sub {
#        my $job = shift;
#        my $opt = shift;
#
#        my $file   = $job;
#        my $output = $file;
#        $output =~ s/fa$/nib/;
#
#        # faToNib - Convert from .fa to .nib format
#        # usage:
#        #   faToNib [options] in.fa out.nib
#        my $cmd = "$kent_bin/faToNib -softMask" . " $file" . " $output";
#        exec_cmd($cmd);
#
#        return;
#    };
#
#    my @jobs = sort @files;
#
#    my $start_time = time;
#    print "\n", "=" x 30, "\n";
#    print "Processing...\n";
#
#    my $run = AlignDB::Run->new(
#        parallel => $parallel,
#        jobs     => \@jobs,
#        code     => $worker,
#    );
#    $run->run;
#
#    print "\n";
#    print "Runtime ", duration( time - $start_time ), ".\n";
#    print "=" x 30, "\n";
#}

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

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
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
        print "Run axtChain...\n";
        my $cmd
            = "$kent_bin/axtChain -minScore=$minScore -linearGap=$linearGap -psl" 
            . " $file"
            . " $dir_target/chr.2bit"
            . " $dir_query/chr.2bit"
            . " $output";
        exec_cmd($cmd);
        print ".chain file generated.\n\n";

        return;
    };

    my @jobs = sort @files;

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

#----------------------------------------------------------#
# chain-net section
#----------------------------------------------------------#
{
    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

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
    $cmd
        = "$kent_bin/chainNet -minSpace=1"
        . " $dir_lav/all.pre.chain"
        . " $dir_target/chr.sizes"
        . " $dir_query/chr.sizes"
        . " $dir_lav/target.chainnet"
        . " /dev/null";
    exec_cmd($cmd);

    # netSyntenic - Add synteny info to net.
    # usage:
    #   netSyntenic in.net out.net
    $cmd
        = "$kent_bin/netSyntenic"
        . " $dir_lav/target.chainnet"
        . " $dir_lav/target.net";
    exec_cmd($cmd);

    # netSplit - Split a genome net file into chromosome net files
    # usage:
    #   netSplit in.net outDir
    $cmd = "$kent_bin/netSplit" . " $dir_lav/target.net" . " $dir_per_chr";
    exec_cmd($cmd);

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

#----------------------------------------------------------#
# netToAxt section
#----------------------------------------------------------#
{
    my @files = File::Find::Rule->file->name('*.net')->in($dir_per_chr);
    printf "\n----%4s .net files to be converted ----\n", scalar @files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file   = $job;
        my $output = $file;
        $output =~ s/net$/axt/;

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
            . " stdout | $kent_bin/axtSort stdin"
            . " $output";
        exec_cmd($cmd);
        print ".axt file generated.\n\n";

        return;
    };

    my @jobs = sort @files;

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

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
    
    >perl part_seq.pl -in t\RM11 -out t\parted -chunk 500000
    >dos2unix *.fa
    >perl bz.pl -dt t\S288C\chr01.fa -dq t\parted -s set01 -dl test -qp
    
    >perl part_seq.pl -in t\S288C -out t\S288C_parted -chunk 500000
    >perl bz.pl -dt t\S288C_parted -dq t\RM11\rm11.fa -s set01 -dl test -tp -p 4

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Blastz will take the first sequence in target fasta file and all sequences in
query fasta file.

So if there are mutiple query files, the program uses the largest one. And all
target files in dir_target will be processed.

So, if there are combined fasta files and multi fasta files coexisting in the
target directory, just delete the axt file matched with the combined fasta
filename.

=cut
