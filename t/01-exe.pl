#!/usr/bin/perl
use strict;
use warnings;

use FindBin;

my $cmd;

$cmd .= "$FindBin::Bin/../blastz";
$cmd .= " $FindBin::Bin/S288C/chr01.fa";
$cmd .= " $FindBin::Bin/RM11/supercontig1_1.fa";
$cmd .= " > chr01.new.lav";

print $cmd, "\n";
system($cmd);