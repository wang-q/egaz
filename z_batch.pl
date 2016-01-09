#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Time::Duration;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

z_batch.pl - bz.pl, lpcna.pl and amp.pl

=head1 SYNOPSIS

    perl z_batch.pl -dt t/S288C -dq t/RM11 -dw . -p 4 -r 1

    perl z_batch.pl -dt t/S288C -dq t/RM11 -dw . -p 4 -r 2-4

=cut

GetOptions(
    'help|?'          => sub { HelpMessage(0) },
    'dir_target|dt=s' => \my $dir_target,
    'dir_query|dq=s'  => \my $dir_query,
    'dir_working|dw=s' => \( my $dir_working = '.' ),
    'parallel|p=i'     => \( my $parallel    = 1 ),
    'run|r=s'          => \( my $task        = 1 ),
    'clean'            => \( my $clean ),
) or HelpMessage(1);

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
    my $t_base = path($dir_target)->basename;
    $t_base =~ s/\..+?$//;
    my $q_base = path($dir_query)->basename;
    $q_base =~ s/\..+?$//;
    my $tq = "${t_base}vs${q_base}";
    $ldir = path( $dir_working, $tq );
    $ldir = $ldir->absolute->stringify;
    if ( $clean and path($ldir)->exists ) {
        path($ldir)->remove_tree;
    }
}

#----------------------------------------------------------#
# dispatch table
#----------------------------------------------------------#
my $dispatch = {
    1 => "perl $FindBin::Bin/bz.pl"
        . " -s set01"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    2 => "perl $FindBin::Bin/bz.pl"
        . " -s set01 --noaxt"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    3 => "perl $FindBin::Bin/lpcna.pl"
        . " -dt $dir_target"
        . " -dq $dir_query"
        . " -dl $ldir"
        . " --parallel $parallel",
    4 => "perl $FindBin::Bin/amp.pl" . " -syn"
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

    if ( $step == 1 or $step == 2 ) {
        path($ldir)->remove_tree;
        path($ldir)->mkpath;
    }

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

__END__
