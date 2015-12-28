#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Graph;
use List::MoreUtils qw(uniq);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(string_to_set);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

merge_node.pl - merge overlapped nodes of paralog graph
    
=head1 SYNOPSIS

    perl merge_node.pl -f <file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     tsv link files
        --output        -o  STR     output   
        --coverage      -c  FLOAT   When larger than this ratio, merge nodes, default is [0.9]       
        --verbose       -v          verbose mode

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'files|f=s'  => \my @files,
    'output|o=s' => \my $output,
    'coverage|c=f' => \( my $coverage = 0.9 ),
    'v|verbose|v' => \my $verbose,
) or HelpMessage(1);

die "Need --file\n" if @files == 0;
for my $file (@files) {
    if ( !path($file)->is_file ) {
        die "--file [$file] doesn't exist\n";
    }
}

if ( !$output ) {
    $output = path( $files[0] )->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.merge.yml";
}

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Paralog graph");

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

#----------------------------#
# Read
#----------------------------#
my $g = Graph->new( directed => 0 );
my %chrs;

# nodes are in the form of "I(+):50-100"
for my $file (@files) {
    $stopwatch->block_message("Loading [$file]");
    my @lines = path($file)->lines( { chomp => 1 } );

    for my $line (@lines) {
        my @nodes = ( split /\t/, $line )[ 0, 1 ];
        for my $node (@nodes) {
            if ( !$g->has_vertex($node) ) {
                $g->add_vertex($node);
                my ( $chr, $set, $strand ) = string_to_set($node);
                $g->set_vertex_attribute( $node, "chr",    $chr );
                $g->set_vertex_attribute( $node, "set",    $set );
                $g->set_vertex_attribute( $node, "strand", $strand );

                $chrs{$chr}++;
                print "Add node $node\n";
            }
        }
    }
    $stopwatch->block_message("Finish loading [$file]");
}

#----------------------------#
# Merge
#----------------------------#
$stopwatch->block_message("Merge nodes");

for my $chr ( sort keys %chrs ) {
    print "Merge nodes in chromosome [$chr]\n";
    my @nodes = sort grep { $g->get_vertex_attribute( $_, "chr" ) eq $chr } $g->vertices;

    for my $i ( 0 .. $#nodes ) {
        my $node_i = $nodes[$i];
        print " " x 4, "Node $i / @{[$#nodes]}\t$node_i\n";
        my $set_i = $g->get_vertex_attribute( $node_i, "set" );
        for my $j ( $i + 1 .. $#nodes ) {
            my $node_j = $nodes[$j];
            my $set_j = $g->get_vertex_attribute( $node_j, "set" );

            my $i_set = $set_i->intersect($set_j);
            if ( $i_set->is_not_empty ) {
                my $coverage_i = $i_set->size / $set_i->size;
                my $coverage_j = $i_set->size / $set_j->size;
                if (    $coverage_i >= $coverage
                    and $coverage_j >= $coverage )
                {
                    $g->add_edge( $nodes[$i], $nodes[$j] );
                    print " " x 8, "Merge with Node $j / @{[$#nodes]}\t$node_j\n";
                }
            }
        }
    }
}

#----------------------------#
# Hash of merge
#----------------------------#
$stopwatch->block_message("Output merged");
my $merged_of = {};
{
    my @cc = $g->connected_components;

    # filter single nodes
    @cc = grep { scalar @{$_} > 1 } @cc;

    for my $c (@cc) {
        print "\n";
        printf "    To merge %s nodes\n", scalar @{$c};
        my $chr = $g->get_vertex_attribute( $c->[0], "chr" );
        my $merge_set = AlignDB::IntSpan->new;
        my @strands;
        my ( $strand, $change );
        for my $node ( @{$c} ) {
            my $set = $g->get_vertex_attribute( $node, "set" );
            $merge_set->add($set);
            push @strands, $g->get_vertex_attribute( $node, "strand" );
        }
        @strands = uniq(@strands);
        if ( @strands == 1 ) {
            print " " x 4, "All nodes have the same strand\n";
            $strand = $strands[0];
            $change = 0;
        }
        else {
            print " " x 4, "Nodes have different strands\n";
            $strand = "+";
            $change = 1;
        }
        my $merge_node = "$chr($strand):" . $merge_set->runlist;

        for my $node ( @{$c} ) {
            my $node_change = 0;
            if ($change) {
                my $node_strand = $g->get_vertex_attribute( $node, "strand" );
                if ( $node_strand ne $strand ) {
                    $node_change = 1;
                }
            }
            $merged_of->{$node} = { node => $merge_node, change => $node_change };
            print " " x 8 . "$node => $merge_node\n";
        }
    }

    DumpFile( $output, $merged_of );
}

$stopwatch->end_message;
exit;

__END__
