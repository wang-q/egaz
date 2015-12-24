#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use Graph;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

merge_node.pl - merge overlapped nodes of paralog graph
    
=head1 SYNOPSIS

    perl merge_node.pl -f <file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     node pair tsv files
        --output        -o  STR     Graph yml dump
        --graph         -g  STR     precomputed graph
        --merge         -m  STR     merged nodes, hashref
        --nonself       -n          skip self match, even for palindrome       
        --verbose       -v          verbose mode
=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'files|f=s'  => \my @files,
    'output|o=s' => \my $output,
    'graph|g=s'  => \my $graph_file,
    'merge|m=s'  => \my $merge_file,
    'nonself|n'  => \my $nonself,
    'verbose|v'  => \my $verbose,
) or HelpMessage(1);

if ( !$output ) {
    $output = basename( $graph_file ? $graph_file : $files[0] );
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.graph.yml";
}

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Paralog graph");

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

my $merged_of;
if ($merge_file) {
    $merged_of = LoadFile($merge_file);
}

# nodes are in "chr1:50-100" or "chr1(+):50-100" form, and with an attribute of intspan object
for my $file (@files) {
    open my $in_fh, "<", $file;
LINE: while ( my $line = <$in_fh> ) {
        print "\n";
        chomp $line;

        my (@nodes) = ( split /\t/, $line )[ 0, 1 ];
        my ($hit_strand) = ( split /\t/, $line )[2];

        # convert to merged node
        for my $node (@nodes) {
            if ( exists $merged_of->{$node} ) {
                printf "%s => %s\n", $node, $merged_of->{$node}{node}
                    if $verbose;
                $node = $merged_of->{$node}{node};
                if ( $merged_of->{$node}{change} ) {
                    $hit_strand = change_strand($hit_strand);
                }
            }
        }

        if ($nonself) {
            if ( $nodes[0] eq $nodes[1] ) {
                print " " x 4, "Same node, next\n" if $verbose;
                next LINE;
            }
        }

        # add node
        for my $node (@nodes) {
            if ( !$g->has_vertex($node) ) {
                $g->add_vertex($node);
                print "Add node $node\n";
            }
        }

        # add edge
        if ( !$g->has_edge(@nodes) ) {
            $g->add_edge(@nodes);
            $g->set_edge_attribute( @nodes, "strand", $hit_strand );

            print join "\t", @nodes, "\n" if $verbose;
            printf "Nodes %d \t Edges %d\n", scalar $g->vertices, scalar $g->edges
                if $verbose;
        }
        else {
            print " " x 4, "Edge exists, next\n" if $verbose;
        }
    }
    $stopwatch->block_message("Finish processing [$file]");
}

DumpFile( $output, $g );

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub string_to_set {
    my $node = shift;

    my ( $chr, $runlist ) = split /:/, $node;
    my $strand = "+";
    if ( $chr =~ /\((.+)\)/ ) {
        $strand = $1;
        $chr =~ s/\(.+\)//;
    }
    my $set = AlignDB::IntSpan->new($runlist);

    return ( $chr, $set, $strand );
}

sub change_strand {
    my $strand = shift;

    if ( $strand eq '+' ) {
        return '-';
    }
    elsif ( $strand eq '-' ) {
        return '+';
    }
    else {
        return $strand;
    }
}

#sub string_to_hashobj {
#    my $node = shift;
#
#    my $set_of = {};
#    my @chr_runlists = grep {defined} split /;/, $node;
#    for my $chr_runlist (@chr_runlists) {
#        my ( $chr_name, $runlist ) = split /:/, $chr_runlist;
#        if ( !exists $set_of->{$chr_name} ) {
#            $set_of->{$chr_name} = AlignDB::IntSpan->new;
#        }
#        $set_of->{$chr_name}->add($runlist);
#    }
#
#    return $set_of;
#}

#sub hashobj_to_string {
#    my $set_of = shift;
#
#    my $string;
#    for my $chr_name ( sort keys %{$set_of} ) {
#        my $runlist = $set_of->{$chr_name}->runlist;
#        $string .= "$chr_name:$runlist;";
#    }
#    $string =~ s/;$//;
#
#    return $string;
#}
#
#sub hashstring_to_hashobj {
#    my $set_of = shift;
#
#    for my $key ( keys %{$set_of} ) {
#        $set_of->{$key} = AlignDB::IntSpan->new( $set_of->{$key} );
#    }
#
#    return;
#}
#
#sub hashobj_to_hashstring {
#    my $set_of = shift;
#
#    for my $key ( keys %{$set_of} ) {
#        $set_of->{$key} = $set_of->{$key}->runlist;
#    }
#
#    return;
#}
#
#sub merge_runlist {
#    my $node0 = shift;
#    my $node1 = shift;
#
#    # transforming to hashobj
#    $node0 = string_to_hashobj($node0);
#    $node1 = string_to_hashobj($node1);
#
#    my %seen;
#    for my $node ( $node0, $node1 ) {
#        for my $key ( sort keys %{$node} ) {
#            $seen{$key}++;
#        }
#    }
#
#    my $new_node = {};
#    for my $key ( sort keys %seen ) {
#        if ( $seen{$key} == 1 ) {
#            $new_node->{$key}
#                = AlignDB::IntSpan->new( $node0->{$key}, $node1->{$key} );
#        }
#        else {
#            $new_node->{$key} = $node0->{$key}->union( $node1->{$key} );
#        }
#    }
#
#    return hashobj_to_string($new_node);
#}
#
#sub string_chr {
#    my $chr_runlist = shift;
#    my ( $chr_name, $runlist ) = split /:/, $chr_runlist;
#
#    return $chr_name;
#}

__END__
