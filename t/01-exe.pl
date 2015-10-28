#!/usr/bin/perl
use strict;
use warnings;

use FindBin;

my $cmd;

$cmd .= "lastz";
$cmd .= " $FindBin::Bin/S288C/chrI.fa";
$cmd .= " $FindBin::Bin/RM11/RM11.fa";
$cmd .= " > chrI.new.lav";

print $cmd, "\n";
system($cmd);
