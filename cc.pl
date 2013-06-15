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

    my @edges = $g->edges_at($node0);
    $g->delete_vertex($node0);

    for my $edge (@edges) {
        my @nodes = @{$edge};
        @nodes = map { $_ eq $node0 ? $node1 : $_ } @nodes;
        $g->add_vertex($node1);
        $g->set_vertex_attribute( $node1, "set", string_to_hashobj($node1) );
        $g->add_edge(@nodes);
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


