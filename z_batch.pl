#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Spec;
use Time::Duration;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# executable file location
my $kent_bin = "~/bin/x86_64";

# dirs
my ( $dir_target, $dir_query, );
my $dir_working = ".";

# running tasks
my $task = "1";

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'           => \$help,
    'man'              => \$man,
    'bin|kent_bin=s'   => \$kent_bin,
    'dt|dir_target=s'  => \$dir_target,
    'dq|dir_query=s'   => \$dir_query,
    'dw|dir_working=s' => \$dir_working,
    'parallel=i'       => \$parallel,
    'r|run=s'          => \$task,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my @tasks;
{
    $task =~ s/\"\'//s;
    if ( AlignDB::IntSpan->valid($task) ) {
        my $set = AlignDB::IntSpan->new($task);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $task;
    }
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $ldir;
{
    my $t_base = basename($dir_target);
    $t_base =~ s/\..+?$//;
    my $q_base = basename($dir_query);
    $q_base =~ s/\..+?$//;
    my $tq = "${t_base}vs${q_base}";
    $ldir = File::Spec->catdir( $dir_working, $tq );
    $ldir = File::Spec->rel2abs($ldir);
}

#----------------------------------------------------------#
# dispatch table
#----------------------------------------------------------#
my $dispatch = {
    1 => "perl $FindBin::Bin/bz.pl"
        . " -s set01 -pb lastz --lastz"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    2 => "perl $FindBin::Bin/bz.pl"
        . " -s set01 --noaxt -pb lastz --lastz"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    3 => "perl $FindBin::Bin/lpcna.pl"
        . " -bin $kent_bin"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    4 => "perl $FindBin::Bin/amp.pl" . " -syn"
        . " -bin $kent_bin"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
};

#----------------------------#
# Run
#----------------------------#
# use the dispatch template to generate $cmd
for my $step (@tasks) {
    my $cmd = $dispatch->{$step};

    next unless $cmd;

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing Step $step\n";
    exec_cmd($cmd);
    print "Finish Step $step\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}
print "\n";

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

perl z_batch.pl -dt t/S288C -dq t/RM11/rm11.fa -dw . -p 4 -r 1

perl z_batch.pl -dt t/S288C -dq t/RM11/rm11.fa -dw . -p 4 -r 2-4
