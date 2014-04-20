#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use File::Spec;
use Path::Class;
use Graph;
use List::Util qw(sum0);
use Set::Scalar;

use FindBin;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $file;
my $size_file;

my $output;
my $outdir;

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

    ($output) = grep {defined} split /\./, $output;
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
print "Load graph $file\n";
my $g = LoadFile($file);
my $gnew = Graph->new( directed => 0 );

#----------------------------#
# create cc
#----------------------------#
print "Get sorted cc\n";
my @cc = cc_sorted($g);

#----------------------------#
# Change strands of nodes based on first node in cc
#----------------------------#
print "Change strands\n";
my $count_old_of = {};
for my $c (@cc) {
    my @nodes = @{$c};
    my $copy  = scalar @nodes;
    $count_old_of->{$copy}++;

    next if $copy == 1;

    print " " x 4, "Copy number of this cc is $copy\n";

    # transform string to array_ref
    my @nodes_ref = map { [ ( string_to_set($_) )[ 0, 1 ], undef ] } @nodes;
    $nodes_ref[0]->[2] = "+";    # set first node to positive strand

    my $assigned  = AlignDB::IntSpan->new(0);
    my $unhandled = AlignDB::IntSpan->new;
    $unhandled->add_pair( 0, $copy - 1 );
    $unhandled->remove($assigned);

    my @edges;                   # store edge pairs. not all nodes connected
    while ( $assigned->size < $copy ) {
        for my $i ( $assigned->elements ) {
            for my $j ( $unhandled->elements ) {
                next if !$g->has_edge( $nodes[$i], $nodes[$j] );
                next
                    if !$g->has_edge_attribute( $nodes[$i], $nodes[$j],
                    "strand" );

                my $edge_strand
                    = $g->get_edge_attribute( $nodes[$i], $nodes[$j],
                    "strand" );
                if ( $edge_strand eq "-" ) {
                    printf " " x 8 . "change strand of %s\n", $nodes[$j];
                    $nodes_ref[$j]->[2] = change_strand( $nodes_ref[$i]->[2] );
                }
                else {
                    $nodes_ref[$j]->[2] = $nodes_ref[$i]->[2];
                }
                $unhandled->remove($j);
                $assigned->add($j);
                push @edges, [ $i, $j ];
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
print "Get sorted new cc\n";
my @cc_new   = cc_sorted($gnew);
my $count_of = {};
for my $c (@cc_new) {
    my $copy = scalar @{$c};
    $count_of->{$copy}++;
}

#----------------------------#
# runlist
#----------------------------#
my @copies = sort { $a <=> $b } keys %{$count_of};
my $set_of = {};
for my $c (@cc) {
    my $copy = scalar @{$c};
    if ( !exists $set_of->{$copy} ) {
        $set_of->{$copy} = {};
    }
    for ( @{$c} ) {
        my ( $chr, $set, $strand ) = string_to_set($_);
        if ( !exists $set_of->{$copy}{$chr} ) {
            $set_of->{$copy}{$chr} = AlignDB::IntSpan->new;
        }
        $set_of->{$copy}{$chr}->add($set);
    }
}

for my $key_i ( keys %{$set_of} ) {
    for my $key_j ( keys %{ $set_of->{$key_i} } ) {
        $set_of->{$key_i}{$key_j} = $set_of->{$key_i}{$key_j}->runlist;
    }
}

DumpFile(
    "$output.cc.yml",
    {   count   => $count_of,
        cc      => \@cc_new,
        runlist => $set_of,
    }
);

DumpFile( "$output.new.graph.yml", $gnew );

$stopwatch->end_message;
exit;

sub cc_sorted {
    my $g  = shift;
    my @cc = $g->connected_components;

    # sort by chromosome order within cc
    for my $c (@cc) {
        my @unsorted = @{$c};

        # start point on chromosomes
        @unsorted = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @unsorted;

        # chromosome name
        @unsorted = map { $_->[0] }
            sort { $a->[1] cmp $b->[1] }
            map { /([\w.]+)\(.\)\:/; [ $_, $1 ] } @unsorted;

        $c = [@unsorted];
    }

    # sort by first node's chromosome order between cc
    @cc = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { $_->[0] =~ /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @cc;

    @cc = map { $_->[0] }
        sort { $a->[1] cmp $b->[1] }
        map { $_->[0] =~ /([\w.]+)\(.\)\:/; [ $_, $1 ] } @cc;

    # sort by nodes number between cc
    @cc = sort { scalar @{$b} <=> scalar @{$a} } @cc;

    return @cc;
}

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

sub set_to_string {
    my $node = shift;

    my $string = $node->[0] . "(" . $node->[2] . "):" . $node->[1]->runlist;

    return $string;
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


