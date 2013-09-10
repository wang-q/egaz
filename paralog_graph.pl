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
my @files;         # node pair tsv files
my $graph_file;    # precomputed graph

my $output;        # Graph yml dump

my $merge_file;    # merged nodes, hashref

my $verbose;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|files=s'  => \@files,
    'g|graph=s'  => \$graph_file,
    'o|output=s' => \$output,
    'm|merge=s'  => \$merge_file,
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
    $output = basename( $graph_file ? $graph_file : $files[0] );
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

my $merged_of;
if ($merge_file) {
    $merged_of = LoadFile($merge_file);
}

# nodes are in "chr1:50-100" form, and with an attribute of intspan object
for my $file (@files) {
    open my $in_fh, "<", $file;
    while ( my $line = <$in_fh> ) {
        chomp $line;

        my @nodes = ( split /\t/, $line )[ 0, 1 ];
        for my $node (@nodes) {

            # convert to merged node
            if ( exists $merged_of->{$node} ) {
                print "$node => @{[$merged_of->{$node}]}\n" if $verbose;
                $node = $merged_of->{$node};
            }

            # add node
            if ( !$g->has_vertex($node) ) {
                $g->add_vertex($node);
                $g->set_vertex_attribute( $node, "set",
                    string_to_hashobj($node) );
                print "Add node $node\n";
            }
        }

        # add edge
        $g->add_edge(@nodes);

        print join "\t", @nodes, "\n" if $verbose;
        print "Nodes @{[scalar $g->vertices]} \t Edges @{[scalar $g->edges]}\n\n"
            if $verbose;
    }
    $stopwatch->block_message("Finish processing [$file]");
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


