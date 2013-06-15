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
my $graph_file;    # precomputed graph

my $output;

my $merge;
my $ratio_cover = 0.9;    # When larger than this ratio, merge runlists

my $verbose;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$file,
    'g|graph=s'  => \$graph_file,
    'o|output=s' => \$output,
    'm|merge'    => \$merge,
    'r|ratio=s'  => \$ratio_cover,
    'v|verbose'  => \$verbose,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Paralog graph");

if ( !$output ) {
    $output = basename($file);
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.graph.yml";
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $g;
if ($graph_file) {
    $g = LoadFile($graph_file);
}
else {
    $g = Graph->new( directed => 0 );
}

# nodes are in "chr1:50-100" form, and with an attribute of intspan object
if ($file) {
    open my $in_fh, "<", $file;
    while ( my $line = <$in_fh> ) {
        chomp $line;

        my @nodes;
        ( $nodes[0], $nodes[1], ) = split /\t/, $line;

        $g->add_vertex( $nodes[0] );
        $g->set_vertex_attribute( $nodes[0], "set",
            string_to_hashobj( $nodes[0] ) );
        $g->add_vertex( $nodes[1] );
        $g->set_vertex_attribute( $nodes[1], "set",
            string_to_hashobj( $nodes[1] ) );
        $g->add_edge( $nodes[0], $nodes[1] );

        print Dump [ "node pair", \@nodes ] if $verbose;
        print "Nodes @{[scalar $g->vertices]} \t Edges @{[scalar $g->edges]}\n"
            if $verbose;
    }
    $stopwatch->block_message("Finish processing [$file]");
}

if ($merge) {
    $stopwatch->block_message("Merge nodes");

REDO: while (1) {
        my @nodes_start = $g->vertices;
        for my $node0 (@nodes_start) {
            my $node0_merge = $g->get_vertex_attribute( $node0, "merge" );
            next if $node0_merge;

            # skip scanned nodes
            $g->set_vertex_attribute( $node0, "merge", 1 );

            my @nodes_merge = find_merge( $g, $node0, $ratio_cover );
            next unless scalar @nodes_merge > 0;
            print Dump [
                "Merge @{[scalar(@nodes_merge) + 1]} runlists",
                [ $node0, @nodes_merge ]
                ]
                if $verbose;

            my $new_node = $node0;
            for my $node1 (@nodes_merge) {
                $new_node = merge_runlist( $new_node, $node1 );
            }

            for my $node2 ( $node0, @nodes_merge ) {
                replace_node( $g, $node2, $new_node, { merge => 1 } );
            }

            print
                "Nodes @{[scalar $g->vertices]} \t Edges @{[scalar $g->edges]}\n"
                if $verbose;
            next REDO;
        }

        if ( scalar $g->vertices == scalar @nodes_start ) {
            last REDO;
        }
    }
}

DumpFile( $output, $g );

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub string_to_hashobj {
    my $node = shift;

    my $set_of = {};
    my @chr_runlists = grep {defined} split /;/, $node;
    for my $chr_runlist (@chr_runlists) {
        my ( $chr_name, $runlist ) = split /:/, $chr_runlist;
        if ( !exists $set_of->{$chr_name} ) {
            $set_of->{$chr_name} = AlignDB::IntSpan->new;
        }
        $set_of->{$chr_name}->add($runlist);
    }

    return $set_of;
}

sub hashobj_to_string {
    my $set_of = shift;

    my $string;
    for my $chr_name ( sort keys %{$set_of} ) {
        my $runlist = $set_of->{$chr_name}->runlist;
        $string .= "$chr_name:$runlist;";
    }
    $string =~ s/;$//;

    return $string;
}

sub hashstring_to_hashobj {
    my $set_of = shift;

    for my $key ( keys %{$set_of} ) {
        $set_of->{$key} = AlignDB::IntSpan->new( $set_of->{$key} );
    }

    return;
}

sub hashobj_to_hashstring {
    my $set_of = shift;

    for my $key ( keys %{$set_of} ) {
        $set_of->{$key} = $set_of->{$key}->runlist;
    }

    return;
}

sub merge_runlist {
    my $node0 = shift;
    my $node1 = shift;

    # transforming to hashobj
    $node0 = string_to_hashobj($node0);
    $node1 = string_to_hashobj($node1);

    my %seen;
    for my $node ( $node0, $node1 ) {
        for my $key ( sort keys %{$node} ) {
            $seen{$key}++;
        }
    }

    my $new_node = {};
    for my $key ( sort keys %seen ) {
        if ( $seen{$key} == 1 ) {
            $new_node->{$key}
                = AlignDB::IntSpan->new( $node0->{$key}, $node1->{$key} );
        }
        else {
            $new_node->{$key} = $node0->{$key}->union( $node1->{$key} );
        }
    }

    return hashobj_to_string($new_node);
}

sub string_chr {
    my $chr_runlist = shift;
    my ( $chr_name, $runlist ) = split /:/, $chr_runlist;

    return $chr_name;
}

sub find_merge {
    my $g     = shift;
    my $node0 = shift;
    my $ratio = shift;

    return unless $g->has_vertex($node0);    # node might be replaced

    my $chr_name  = string_chr($node0);
    my $node0_set = string_to_hashobj($node0);

    my @nodes;
    for my $node1 ( $g->vertices ) {
        my $node1_merge = $g->get_vertex_attribute( $node1, "merge" );
        next if $node1_merge;                # skip merged nodes

        my $node1_set = $g->get_vertex_attribute( $node1, "set" );
        next unless exists $node1_set->{$chr_name};   # skip nodes on other chrs

        next if $node1 eq $node0;                     # skip the same node

        my $i_set
            = $node0_set->{$chr_name}->intersect( $node1_set->{$chr_name} );
        if ( $i_set->is_not_empty ) {
            my $node0_coverage = $i_set->size / $node0_set->{$chr_name}->size;
            my $node1_coverage = $i_set->size / $node1_set->{$chr_name}->size;
            if ( $node0_coverage >= $ratio and $node1_coverage >= $ratio ) {
                push @nodes, $node1;
            }
        }
    }

    return @nodes;
}

sub find_intersect_with_graph {
    my $g     = shift;
    my $node0 = shift;
    my $ratio = shift;

    return $node0 if $g->has_vertex($node0);

    my $node0_set = string_to_hashobj($node0);

    my @nodes = $g->vertices;

    for my $node1 (@nodes) {
        my $node1_set = $g->get_vertex_attribute( $node1, "set" );

        my %seen;
        for my $node ( $node0_set, $node1_set ) {
            for my $key ( sort keys %{$node} ) {
                $seen{$key}++;
            }
        }

        for my $key ( sort grep { $seen{$_} > 1 } keys %seen ) {
            my $i_set = $node0_set->{$key}->intersect( $node1_set->{$key} );
            if ( $i_set->is_not_empty ) {
                my $node0_coverage = $i_set->size / $node0_set->{$key}->size;
                my $node1_coverage = $i_set->size / $node1_set->{$key}->size;
                if ( $node0_coverage >= $ratio and $node1_coverage >= $ratio ) {
                    return $node1;
                }
            }
        }
    }

    return;
}

sub replace_node {
    my $g     = shift;
    my $node0 = shift;
    my $node1 = shift;
    my $attr  = shift;    # additional attr for $node1

    my @edges = $g->edges_at($node0);
    $g->delete_vertex($node0);

    for my $edge (@edges) {
        my @nodes = @{$edge};
        @nodes = map { $_ eq $node0 ? $node1 : $_ } @nodes;
        $g->add_vertex($node1);
        $g->add_edge(@nodes);
        $g->set_vertex_attribute( $node1, "set", string_to_hashobj($node1) );
        if ( $attr and ref $attr eq "HASH" ) {
            for my $key ( keys %{$attr} ) {
                $g->set_vertex_attribute( $node1, $key, $attr->{$key} );
            }
        }
    }

    return;
}

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


