#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use IO::Zlib;
use Graph;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $file;

my $output;

my $verbose;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$file,
    'o|output=s' => \$output,
    'v|verbose'  => \$verbose,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$file]");

if ( !$output ) {
    $output = basename($file);

    #($output) = grep {defined} split /\./, $output;
    $output = "$output.cc";
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $g = LoadFile($file);

my @cc = $g->connected_components;

# sort by chromosome order within cc
for my $c (@cc) {
    $c = [ sort { $a cmp $b } @{$c} ];
}

# sort by first node's chromosome order between cc
@cc = sort { $a->[0] cmp $b->[0] } @cc;

# sort by nodes number between cc
@cc = sort { scalar @{$b} <=> scalar @{$a} } @cc;

# count
my $count = {};
for my $c (@cc) {
    my $size = scalar @{$c};
    $count->{$size}++;
}
DumpFile( "$output.yml", { count => $count, cc => \@cc, } );

# write circos link file
open my $link2_fh, ">", "$output.link2.txt";
open my $linkN_fh, ">", "$output.linkN.txt";
for my $c (@cc) {
    my $size = scalar @{$c};
    next if $size < 2;
    for my $idx1 ( 0 .. $size - 1 ) {
        for my $idx2 ( $idx1 + 1 .. $size - 1 ) {
            my @fields;
            for ( $idx1, $idx2 ) {

                # chr1:244401-246892
                # ["chr1", 244_401, 246_892]
                push @fields, split /\:|\-/, $c->[$_];
            }

            if ( $size > 2 ) {
                print {$linkN_fh} join( " ", @fields ), "\n";
            }
            else {
                print {$link2_fh} join( " ", @fields ), "\n";
            }
        }
    }
}
close $link2_fh;
close $linkN_fh;

$stopwatch->end_message;
exit;

__END__


=head1 NAME

    gather_info_axt.pl - 

=head1 SYNOPSIS

    gather_info_axt.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -t, --ft            target file (output)
        -m, --fm            merge file
        -f, --fields        fields

    perl gather_info_axt.pl -f example.match.tsv

=cut


