#!/usr/bin/perl
use strict;
use warnings;

use FindBin;

my $cmd;

$cmd .= "perl $FindBin::Bin/../bz.pl";
$cmd .= " -p 4";
$cmd .= " -dt $FindBin::Bin/S288C";
$cmd .= " -dq $FindBin::Bin/RM11/RM11.fa";
$cmd .= " -dl $FindBin::Bin/S288CvsRM11_p";

print $cmd, "\n";
system($cmd);
