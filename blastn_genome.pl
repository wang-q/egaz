#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use MCE;
use MCE::Flow Sereal => 1;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

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

    perl blastn_genome.pl -f <fasta file> -g <genome file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     blast report file
        --genome        -g  STR     reference genome file
        --coverage      -c  FLOAT   coverage of identical matches, default is [0.9]       
        --output        -o  STR     output
        --parallel      -p  INT     default is [8]
        --chunk_size        INT     default is [500000]

=cut

GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s' => \( my $file ),
    'coverage|c=f' => \( my $coverage   = 0.95 ),
    'genome|g=s'   => \( my $genome ),
    'output|o=s'   => \( my $output ),
    'parallel|p=i' => \( my $parallel   = 8 ),
    'chunk_size=i' => \( my $chunk_size = 500000 ),
) or Getopt::Long::HelpMessage(1);

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
# Parse blast reports
#----------------------------------------------------------#
$stopwatch->block_message("Parse reports");

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $wid = MCE->wid;
    print "* Process task [$chunk_id] by worker #$wid\n";

    my @lines = @{$chunk_ref};
    my %heads;
    for my $line (@lines) {
        next if $line =~ /^#/;
        chomp $line;

        # qseqid sseqid qstart qend sstart send qlen slen nident
        my @fields = grep {defined} split /\s+/, $line;
        if ( @fields != 9 ) {
            print "Fields error: $line\n";
            next;
        }

        my $query_name = $fields[0];
        my $hit_name   = $fields[1];

        my $query_length    = $fields[6];
        my $identical_match = $fields[8];

        my ( $h_start, $h_end ) = ( $fields[4], $fields[5] );
        if ( $h_start > $h_end ) {
            ( $h_start, $h_end ) = ( $h_end, $h_start );
        }

        my $query_coverage = $identical_match / $query_length;

        if ( $query_coverage >= $coverage ) {
            my $head = "$hit_name(+):$h_start-$h_end";
            $heads{$head}++;
        }
    }

    printf " " x 4 . "Gather %d pieces\n", scalar keys %heads;
    printf " " x 4 . "Last query name is %s\n", ( split /\s+/, $lines[-1] )[0];
    MCE->gather(%heads);
};

MCE::Flow::init {
    chunk_size  => $chunk_size,
    max_workers => $parallel,
};
my %locations = mce_flow_f $worker, $file;
MCE::Flow::finish;

$stopwatch->block_message( "Finish blasting", 1 );

{
    #----------------------------#
    # remove locations fully contained by others
    #----------------------------#
    $stopwatch->block_message("Convert nodes to sets");
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

    $stopwatch->block_message("Sort by positions at chromosomes");

    # end positions
    for my $chr ( sort keys %{$nodes_of_chr} ) {
        my @temps = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { /\:(?:\d+)\-(\d+)/; [ $_, $1 ] } @{ $nodes_of_chr->{$chr} };
        $nodes_of_chr->{$chr} = \@temps;
    }

    # start positions
    for my $chr ( sort keys %{$nodes_of_chr} ) {
        my @temps = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { /\:(\d+)/; [ $_, $1 ] } @{ $nodes_of_chr->{$chr} };
        $nodes_of_chr->{$chr} = \@temps;
    }

    $stopwatch->block_message("Overlaps between sets");
    my %to_remove;
    for my $chr ( sort keys %{$nodes_of_chr} ) {
        my @nodes = @{ $nodes_of_chr->{$chr} };

        my $vicinity = 10;
        for my $idx ( 0 .. $#nodes - $vicinity ) {

            for my $i ( 0 .. $vicinity - 1 ) {
                for my $j ( $i .. $vicinity - 1 ) {
                    my $node_i = $nodes[ $idx + $i ];
                    my $set_i  = $set_of{$node_i};

                    my $node_j = $nodes[ $idx + $j ];
                    my $set_j  = $set_of{$node_j};

                    if ( $set_i->larger_than($set_j) ) {
                        $to_remove{$node_j}++;
                    }
                    elsif ( $set_j->larger_than($set_i) ) {
                        $to_remove{$node_i}++;
                    }
                }
            }
        }
    }

    printf " " x 4 . "Remove [%d] nested nodes\n", scalar keys %to_remove;
    my @sorted;
    for my $chr ( sort keys %{$nodes_of_chr} ) {
        for my $node ( @{ $nodes_of_chr->{$chr} } ) {
            if ( !exists $to_remove{$node} ) {
                push @sorted, $node;
            }
        }
    }

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

sub string_to_set {
    my $node = shift;

    my ( $chr_part, $runlist ) = split /:/, $node;

    my $head_qr = qr{
        (?:(?P<name>[\w_]+)\.)?
        (?P<chr_name>[\w-]+)
        (?:\((?P<chr_strand>.+)\))?
    }xi;
    $chr_part =~ $head_qr;

    my $chr    = $+{chr_name};
    my $strand = "+";
    if ( defined $+{chr_strand} ) {
        $strand = $+{chr_strand};
    }
    my $set = AlignDB::IntSpan->new($runlist);

    return ( $chr, $set, $strand );
}

sub get_seq_faidx {
    my $genome   = shift;
    my $location = shift;    # I:1-100

    my $cmd = sprintf "samtools faidx %s %s", $genome, $location;
    open my $fh_pipe, '-|', $cmd;

    my $seq;
    while ( my $line = <$fh_pipe> ) {
        chomp $line;
        if ( $line =~ /^[\w-]+/ ) {
            $seq .= $line;
        }
    }
    close $fh_pipe;

    return $seq;
}

sub exec_cmd {
    my $string = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $string , "\n";
    print "-" x 30, "\n";

    system $string;
}

__END__
