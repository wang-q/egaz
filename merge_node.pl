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

use MCE;
use MCE::Flow Sereal => 1;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(string_to_set);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

merge_node.pl - merge overlapped nodes via overlapping graph
    
=head1 SYNOPSIS

    perl merge_node.pl -f <file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     tsv link files
        --output        -o  STR     output   
        --coverage      -c  FLOAT   When larger than this ratio, merge nodes, default is [0.9]       
        --parallel      -p  INT     default is [8]
        --verbose       -v          verbose mode

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'files|f=s'  => \my @files,
    'output|o=s' => \my $output,
    'coverage|c=f' => \( my $coverage = 0.9 ),
    'parallel|p=i' => \( my $parallel = 8 ),
    'v|verbose|v'  => \my $verbose,
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
# nodes are in the form of "I(+):50-100"
print $stopwatch->block_message("Convert nodes to sets");
my $graph_of_chr = {};
for my $file (@files) {
    $stopwatch->block_message("Loading [$file]");
    my @lines = path($file)->lines( { chomp => 1 } );

    for my $line (@lines) {
        my @nodes = ( split /\t/, $line )[ 0, 1 ];
        for my $node (@nodes) {
            my ( $chr, $set, $strand ) = string_to_set($node);
            if ( !exists $graph_of_chr->{$chr} ) {
                $graph_of_chr->{$chr} = Graph->new( directed => 0 );
            }
            if ( !$graph_of_chr->{$chr}->has_vertex($node) ) {

                $graph_of_chr->{$chr}->add_vertex($node);
                $graph_of_chr->{$chr}->set_vertex_attribute( $node, "chr",    $chr );
                $graph_of_chr->{$chr}->set_vertex_attribute( $node, "set",    $set );
                $graph_of_chr->{$chr}->set_vertex_attribute( $node, "strand", $strand );

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

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $chr = $chunk_ref->[0];
    my $wid = MCE->wid;
    print "* Process chromosome [$chr] by worker #$wid\n";

    my $g     = $graph_of_chr->{$chr};
    my @nodes = sort $g->vertices;
    my @edges;
    for my $i ( 0 .. $#nodes ) {
        my $node_i = $nodes[$i];
        printf " " x 4 . "Node %d / %d\t%s\n", $i, $#nodes, $node_i;
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
                    push @edges, [ $nodes[$i], $nodes[$j] ];
                    printf " " x 8 . "Merge with Node %d / %d\t%s\n", $j, $#nodes, $node_j;
                }
            }
        }
    }
    printf "Gather %d edges\n", scalar @edges;

    MCE->gather( $chr, \@edges );
};
MCE::Flow::init {
    chunk_size  => 1,
    max_workers => $parallel,
};
my %edge_of = mce_flow $worker, ( sort keys %{$graph_of_chr} );
MCE::Flow::finish;

$stopwatch->block_message("Add edges");
for my $chr ( sort keys %{$graph_of_chr} ) {
    print " " x 4 . "Add edges to chromosome [$chr]\n";
    for my $edge ( @{ $edge_of{$chr} } ) {
        $graph_of_chr->{$chr}->add_edge( @{$edge} );
    }
}

#----------------------------#
# Hash of merge
#----------------------------#
$stopwatch->block_message("Output merged");
my $merged_of = {};
for my $chr ( sort keys %{$graph_of_chr} ) {
    my $g  = $graph_of_chr->{$chr};
    my @cc = $g->connected_components;

    # filter single nodes
    @cc = grep { scalar @{$_} > 1 } @cc;

    for my $c (@cc) {
        print "\n";
        printf " " x 4 . "Merge %s nodes\n", scalar @{$c};
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
            print " " x 4 . "All nodes have the same strand\n";
            $strand = $strands[0];
            $change = 0;
        }
        else {
            print " " x 4 . "Nodes have different strands\n";
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
