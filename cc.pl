#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Graph;
use List::MoreUtils qw(minmax);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(string_to_set set_to_string change_strand sort_cc);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

cc.pl - connected components of paralog graph

=head1 SYNOPSIS

    perl cc.pl -f <file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     file
        --ratio         -r  FLOAT   break links between nodes if lengths differences
                                    large than this ratio, default is [0.8]
        --output        -o  STR     output

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \( my $file ),
    'ratio|r=f' => \( my $ratio = 0.8 ),
    'output|o=s' => \( my $output ),
) or HelpMessage(1);

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
}

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$file]");

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$stopwatch->block_message("Load graph [$file]");

my $g = LoadFile($file);
my $gnew = Graph->new( directed => 0 );

#----------------------------#
# create cc
#----------------------------#
$stopwatch->block_message("Get sorted cc");

my @cc = sort_cc($g->connected_components);

#----------------------------#
# Change strands of nodes based on first node in cc
#----------------------------#
$stopwatch->block_message("Change strands and break bad links");
for my $c (@cc) {
    my @nodes = @{$c};
    my $copy  = scalar @nodes;

    next if $copy == 1;

    print "Copy number of this cc is $copy\n";

    # transform string to array_ref
    my @nodes_ref = map { [ ( string_to_set($_) )[ 0, 1 ], undef ] } @nodes;

    # set first node to positive strand
    $nodes_ref[0]->[2] = "+";

    my $assigned  = AlignDB::IntSpan->new(0);
    my $unhandled = AlignDB::IntSpan->new;
    $unhandled->add_pair( 0, $copy - 1 );
    $unhandled->remove($assigned);

    my @edges;    # store edge pairs. not all nodes connected
    while ( $assigned->size < $copy ) {
        for my $i ( $assigned->elements ) {
            for my $j ( $unhandled->elements ) {
                next if !$g->has_edge( $nodes[$i], $nodes[$j] );
                next
                    if !$g->has_edge_attribute( $nodes[$i], $nodes[$j], "strand" );

                my $edge_strand = $g->get_edge_attribute( $nodes[$i], $nodes[$j], "strand" );
                if ( $edge_strand eq "-" ) {
                    printf " " x 4 . "change strand of %s\n", $nodes[$j];
                    $nodes_ref[$j]->[2] = change_strand( $nodes_ref[$i]->[2] );
                }
                else {
                    $nodes_ref[$j]->[2] = $nodes_ref[$i]->[2];
                }
                $unhandled->remove($j);
                $assigned->add($j);

                my $l_i = $nodes_ref[$i]->[1]->size;
                my $l_j = $nodes_ref[$j]->[1]->size;
                my ( $l_min, $l_max ) = minmax( $l_i, $l_j );
                my $diff_ratio = sprintf "%.3f", $l_min / $l_max;

                if ( $diff_ratio < $ratio ) {
                    printf " " x 4 . "Break links between %s %s\n", $nodes[$i], $nodes[$j];
                    printf " " x 4 . "Ratio[%s]\tMin [%s]\tMax[%s]\n", $diff_ratio, $l_min, $l_max;
                }
                else {
                    push @edges, [ $i, $j ];
                }
            }
        }
    }

    # transform array_ref back to string
    my @nodes_new = map { set_to_string($_) } @nodes_ref;

    for my $edge (@edges) {
        $gnew->add_edge( $nodes_new[ $edge->[0] ], $nodes_new[ $edge->[1] ] );
    }
}

#----------------------------#
# recreate cc
#----------------------------#
$stopwatch->block_message("Get sorted new cc");
my @cc_new = sort_cc($gnew->connected_components);

$stopwatch->block_message("Write [$output.cc.raw.yml]");
DumpFile( "$output.cc.raw.yml", \@cc_new );

$stopwatch->end_message;
exit;

__END__
