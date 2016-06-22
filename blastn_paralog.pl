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
use List::Util qw(max);
use List::MoreUtils qw(uniq);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

blastn_paralog.pl - Link paralog sequences
    
=head1 SYNOPSIS

    perl blastn_paralog.pl -f <fasta file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     blast report file
        --coverage      -c  FLOAT   default is [0.9]       
        --output        -o  STR     output
        --parallel      -p  INT     default is [8]
        --chunk_size        INT     default is [500000]

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \( my $file ),
    'coverage|c=f' => \( my $coverage   = 0.9 ),
    'output|o=s'   => \( my $output ),
    'parallel|p=i' => \( my $parallel   = 8 ),
    'chunk_size=i' => \( my $chunk_size = 500000 ),
) or HelpMessage(1);

if ( !defined $file ) {
    die "Need --file\n";
}
elsif ( !path($file)->is_file ) {
    die "--file [$file] doesn't exist\n";
}

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.blast.tsv";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Link paralogs...");

#----------------------------------------------------------#
# Blast
#----------------------------------------------------------#
$stopwatch->block_message("Parse reports");

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $wid = MCE->wid;
    print "* Process task [$chunk_id] by worker #$wid\n";

    my @lines = @{$chunk_ref};
    my @links;
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
        next if $query_name eq $hit_name;

        my $query_length = $fields[6];
        my $hit_length   = $fields[7];
        my $max_length   = max( $query_length, $hit_length );
        next if $query_length / $max_length < $coverage;
        next if $hit_length / $max_length < $coverage;

        my $identical_match = $fields[8];
        next if $identical_match / $max_length < $coverage;

        my ( $h_start, $h_end ) = ( $fields[4], $fields[5] );
        my $strand = "+";
        if ( $h_start > $h_end ) {
            ( $h_start, $h_end ) = ( $h_end, $h_start );
            $strand = "-";
        }

        my $link = join "\t", $query_name, $hit_name, $strand;
        push @links, $link;
    }

    printf "Gather %d links\n", scalar @links;
    MCE->gather(@links);
};

MCE::Flow::init {
    chunk_size  => $chunk_size,
    max_workers => $parallel,
};
my @all_links = mce_flow_f $worker, $file;
MCE::Flow::finish;

#----------------------------------------------------------#
# Write
#----------------------------------------------------------#
$stopwatch->block_message("Remove duplicated links");
@all_links = uniq(@all_links);

path($output)->spew( map {"$_\n"} @all_links );

$stopwatch->end_message;

exit;

__END__
