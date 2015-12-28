#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Bio::SearchIO;
use Set::Scalar;

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

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $file,
    'view|m=s'     => \( my $alignment_view = 0 ),
    'identity|i=i' => \( my $identity       = 90 ),
    'coverage|c=f' => \( my $coverage       = 0.95 ),
    'genome|g=s'   => \my $genome,
    'output|o=s'   => \my $output,
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

my $view_name = {
    0 => "blast",         # Pairwise
    7 => "blastxml",      # BLAST XML
    9 => "blasttable",    # Hit Table
};
my $result_format = $view_name->{$alignment_view};

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

my %locations;    # genome locations
open my $blast_fh, '<', $file;
my $searchio = Bio::SearchIO->new(
    -format => $result_format,
    -fh     => $blast_fh,
);

QUERY: while ( my $result = $searchio->next_result ) {
    my $query_name = $result->query_name;

    # blasttable don't have $result->query_length
    my $query_info   = decode_header($query_name);
    my $query_length = $query_info->{chr_end} - $query_info->{chr_start} + 1;
    print "Name $query_name\tLength $query_length\n";
    while ( my $hit = $result->next_hit ) {
        my $hit_name = $hit->name;

        while ( my $hsp = $hit->next_hsp ) {
            my $query_set = AlignDB::IntSpan->new;

            # process the Bio::Search::HSP::HSPI object
            my $hsp_length = $hsp->length('query');

            # use "+" for default strand
            # -1 = Minus strand, +1 = Plus strand
            my ( $query_strand, $hit_strand ) = $hsp->strand("list");
            my $hsp_strand = "+";
            if ( $query_strand + $hit_strand == 0 and $query_strand != 0 ) {
                $hsp_strand = "-";
            }

            my ( $q_start, $q_end ) = $hsp->range('query');
            if ( $q_start > $q_end ) {
                ( $q_start, $q_end ) = ( $q_end, $q_start );
            }
            $query_set->add_range( $q_start, $q_end );

            my ( $h_start, $h_end ) = $hsp->range('hit');
            if ( $h_start > $h_end ) {
                ( $h_start, $h_end ) = ( $h_end, $h_start );
            }

            my $query_coverage = $query_set->size / $query_length;
            my $hsp_identity   = $hsp->percent_identity;

            #print Dump {
            #    hit_name     => $hit_name,
            #    hsp_identity => $hsp_identity,
            #    q_start      => $q_start,
            #    q_end        => $q_end,
            #    h_start      => $h_start,
            #    h_end        => $h_end,
            #    query_strand => $query_strand,
            #    hit_strand   => $hit_strand,
            #    hsp_strand   => $hsp_strand,
            #};

            if ( $query_coverage >= $coverage and $hsp_identity >= $identity ) {
                my $head = "$hit_name(+):$h_start-$h_end";

                if (    $hit_name eq $query_info->{chr_name}
                    and $h_start eq $query_info->{chr_start}
                    and $h_end eq $query_info->{chr_end} )
                {
                    print " " x 4 . "Find itself\n";
                }
                else {

                    print " " x 4 . "Find new [$head]\n";
                    $locations{$head}++;
                }
            }
        }
    }
}

close $blast_fh;

{
    #----------------------------#
    # remove locations fully contained by others
    #----------------------------#
    print $stopwatch->block_message("Merge nested locations");
    my %chrs;
    my %set_of;
    for my $node ( keys %locations ) {
        my ( $chr, $set, $strand ) = string_to_set($node);
        $chrs{$chr}++;
        $set_of{$node} = { chr => $chr, set => $set };
    }

    my $to_remove = Set::Scalar->new;
    for my $chr ( sort keys %chrs ) {
        my @nodes = sort grep { $set_of{$_}->{chr} eq $chr } keys %locations;

        for my $i ( 0 .. $#nodes ) {
            my $node_i = $nodes[$i];
            my $set_i  = $set_of{$node_i}->{set};
            for my $j ( $i + 1 .. $#nodes ) {
                my $node_j = $nodes[$j];
                my $set_j  = $set_of{$node_j}->{set};

                if ( $set_i->larger_than($set_j) ) {
                    $to_remove->insert($node_j);
                }
                elsif ( $set_j->larger_than($set_i) ) {
                    $to_remove->insert($node_i);
                }
            }
        }
    }

    #----------------------------#
    # sort heads
    #----------------------------#
    $stopwatch->block_message("Sort locations");
    my @sorted = map { $to_remove->has($_) ? () : $_ } keys %locations;

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
