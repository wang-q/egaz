#!/usr/bin/perl
use strict;
use warnings;

use FindBin;

my $cmd;

$cmd .= "perl $FindBin::Bin/../bz.pl";
$cmd .= " -dt $FindBin::Bin/S288C";
$cmd .= " -dq $FindBin::Bin/RM11";
$cmd .= " -dl $FindBin::Bin/S288CvsRM11_d";

print $cmd, "\n";
system($cmd);
