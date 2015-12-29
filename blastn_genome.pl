#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use MCE;
use MCE::Flow Sereal => 1;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(string_to_set get_seq_faidx decode_header);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

blastn_genome.pl - Get more paralog pieces by genomic blasting
    
=head1 SYNOPSIS

    perl blastn_genome.pl -f <blast result file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     blast result file
        --view          -M  STR     blast output format, default is [0]
                                    `blastall -m`
                                    0 => "blast",         # Pairwise
                                    7 => "blastxml",      # BLAST XML
                                    9 => "blasttable",    # Hit Table
        --identity      -i  INT     default is [90]
        --coverage      -c  FLOAT   default is [0.95]       
        --genome        -g  STR     reference genome file
        --output        -o  STR     output
        --parallel      -p  INT     default is [8]
        --chunk_size        INT     default is [1000000]

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \( my $file ),
    'identity|i=i' => \( my $identity   = 90 ),
    'coverage|c=f' => \( my $coverage   = 0.95 ),
    'genome|g=s'   => \( my $genome ),
    'output|o=s'   => \( my $output ),
    'parallel|p=i' => \( my $parallel   = 8 ),
    'chunk_size=i' => \( my $chunk_size = 1000000 ),
) or HelpMessage(1);

if ( !defined $file ) {
    die "Need --file\n";
}
elsif ( !path($file)->is_file ) {
    die "--file [$file] doesn't exist\n";
}

if ( !defined $genome ) {
    die "--genome is needed\n";
}
elsif ( !path($genome)->is_file ) {
    die "--genome doesn't exist\n";
}

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.bg.fasta";
}
path($output)->remove;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Find paralogs...");

#----------------------------------------------------------#
# load blast reports
#----------------------------------------------------------#
$stopwatch->block_message("load blast reports");

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $wid = MCE->wid;
    print "* Process task [$chunk_id] by worker #$wid\n";

    my @lines = @{$chunk_ref};
    my %heads;
    for my $line (@lines) {
        next if $line =~ /^#/;
        chomp $line;

        # Query id, Subject id, %identity, alignment length, mismatches, gap openings,
        # q.start, q.end, s. start, s. end, e-value, bit score
        my @fields = grep {defined} split /\s+/, $line;
        if ( @fields != 12 ) {
            print "Fields error: $line\n";
            next;
        }

        my $query_name = $fields[0];
        my $hit_name   = $fields[1];

        my $query_info   = decode_header($query_name);
        my $query_length = $query_info->{chr_end} - $query_info->{chr_start} + 1;

        my ( $q_start, $q_end ) = ( $fields[6], $fields[7] );
        if ( $q_start > $q_end ) {
            ( $q_start, $q_end ) = ( $q_end, $q_start );
        }
        my $query_set = AlignDB::IntSpan->new;
        $query_set->add_pair( $q_start, $q_end );

        my ( $h_start, $h_end ) = ( $fields[8], $fields[9] );
        if ( $h_start > $h_end ) {
            ( $h_start, $h_end ) = ( $h_end, $h_start );
        }

        my $query_coverage = $query_set->size / $query_length;
        my $hsp_identity   = $fields[2];

        #print Dump {
        #    hit_name       => $hit_name,
        #    hsp_identity   => $hsp_identity,
        #    q_start        => $q_start,
        #    q_end          => $q_end,
        #    h_start        => $h_start,
        #    h_end          => $h_end,
        #    query_coverage => $query_coverage,
        #    hsp_identity   => $hsp_identity,
        #};

        if ( $query_coverage >= $coverage and $hsp_identity >= $identity ) {
            my $head = "$hit_name(+):$h_start-$h_end";
            $heads{$head}++;
        }
    }

    printf "Gather %d pieces\n", scalar keys %heads;
    MCE->gather(%heads);
};

MCE::Flow::init {
    chunk_size  => $chunk_size,
    max_workers => $parallel,
};
my %locations = mce_flow_f $worker, $file;
MCE::Flow::finish;

$stopwatch->block_message( "Finish parse blast results", 1 );

{
    #----------------------------#
    # remove locations fully contained by others
    #----------------------------#
    print $stopwatch->block_message("Convert nodes to sets");
    my $nodes_of_chr = {};
    my %set_of;
    for my $node ( keys %locations ) {
        my ( $chr, $set, $strand ) = string_to_set($node);
        if ( !exists $nodes_of_chr->{$chr} ) {
            $nodes_of_chr->{$chr} = [];
        }
        push @{ $nodes_of_chr->{$chr} }, $node;
        $set_of{$node} = $set;
    }

    print $stopwatch->block_message("Overlaps between sets");
    my $worker = sub {
        my ( $self, $chunk_ref, $chunk_id ) = @_;

        my $chr = $chunk_ref->[0];
        my $wid = MCE->wid;
        print "* Process chromosome [$chr] by worker #$wid\n";

        my @nodes = @{ $nodes_of_chr->{$chr} };
        my %seen;
        for my $i ( 0 .. $#nodes ) {
            my $node_i = $nodes[$i];
            my $set_i  = $set_of{$node_i};
            for my $j ( $i + 1 .. $#nodes ) {
                my $node_j = $nodes[$j];
                my $set_j  = $set_of{$node_j};

                if ( $set_i->larger_than($set_j) ) {
                    $seen{$node_j}++;
                }
                elsif ( $set_j->larger_than($set_i) ) {
                    $seen{$node_i}++;
                }
            }
        }
        MCE->gather(%seen);
    };
    MCE::Flow::init {
        chunk_size  => 1,
        max_workers => $parallel,
    };
    my %to_remove = mce_flow $worker, ( sort keys %{$nodes_of_chr} );
    MCE::Flow::finish;

    #----------------------------#
    # sort heads
    #----------------------------#
    $stopwatch->block_message("Sort locations");
    my @sorted = map { exists $to_remove{$_} ? () : $_ } keys %locations;

    # start point on chromosomes
    @sorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @sorted;

    # chromosome name
    @sorted = map { $_->[0] }
        sort { $a->[1] cmp $b->[1] }
        map { /([\w.]+)\(.\)\:/; [ $_, $1 ] } @sorted;

    #----------------------------#
    # write fasta
    #----------------------------#
    $stopwatch->block_message("Write outputs [$output]");
    for my $node (@sorted) {
        my ( $chr, $set, $strand ) = string_to_set($node);
        my $location = "$chr:" . $set->runlist;
        my $seq = get_seq_faidx( $genome, $location );
        path($output)->append(">$node\n");
        path($output)->append("$seq\n");
    }
}

$stopwatch->end_message;

exit;

__END__
