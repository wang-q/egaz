#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Tie::IxHash;
use List::MoreUtils qw(uniq);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(get_seq_faidx decode_header);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

genome_locations.pl - Get exact genome locations of sequence pieces in a fasta file

=head1 SYNOPSIS

    perl blastn_genome_locations.pl -f <blast result file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     query fasta file
        --genome        -g  STR     reference genome file
        --length        -l  INT     Match length of sparsemem, default is [50]
        --coverage      -c  FLOAT   default is [0.9]
        --hole              INT     Fill holes (snp/indel) less than this length, default is [10]
        --output        -o  STR     output
        --debug                     Write mem.yml for debugging

=head1 REQUIREMENTS

    brew tap homebrew/science
    brew tap wang-q/tap
    brew install sparsemem samtools

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \my $file,
    'genome|g=s' => \my $genome,
    'coverage|c=f' => \( my $coverage     = 0.9 ),
    'hole=i'       => \( my $hole_length  = 10 ),
    'length=i'     => \( my $match_length = 50 ),
    'output|o=s'   => \( my $output ),
    'debug'        => \( my $debug ),
) or HelpMessage(1);

if ( !$output ) {
    $output = path($file)->basename( ".fasta", ".fas", ".fa" );
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.gl.fasta";
}
path($output)->remove;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Get locations...");

#----------------------------------------------------------#
# load blast reports
#----------------------------------------------------------#
$stopwatch->block_message("Generate sparsemem reports");

my $mem_result = run_sparsemem( $file, $genome, $match_length );

tie my %info_of, "Tie::IxHash";
my $cur_name;
my %match_of;
for my $line ( split /\n/, $mem_result ) {
    if ( $line =~ /^\>/ ) {
        $line =~ s/^\>\s*//;
        chomp $line;
        $cur_name = $line;
    }
    elsif ( $line =~ /^\s*[\w-]+/ ) {
        chomp $line;

        #   XII    1065146         1      6369
        # ID
        # the position in the reference sequence
        # the position in the query sequence
        # the length of the match
        my @fields = grep {/\S+/} split /\s+/, $line;

        if ( @fields != 4 ) {
            warn "    Line [$line] seems wrong.\n";
            next;
        }
        if ( !exists $info_of{$cur_name} ) {
            $info_of{$cur_name} = {};
        }

        if ( !exists $info_of{$cur_name}->{ $fields[0] } ) {
            $info_of{$cur_name}->{ $fields[0] } = AlignDB::IntSpan->new;
        }

        $info_of{$cur_name}->{ $fields[0] }->add_pair( $fields[1], $fields[1] + $fields[3] - 1 );
    }
    else {    # Blank line, do nothing
    }
}

my @locations;    #  genome locations
my %count_of;     # hit counts
for my $seq_name ( keys %info_of ) {
    my $real_name = $seq_name;
    $real_name =~ s/ Reverse//;
    my $info   = decode_header($real_name);
    my $length = $info->{chr_end} - $info->{chr_start} + 1;

    for my $chr ( keys $info_of{$seq_name} ) {
        if ( $hole_length > 0 ) {
            $info_of{$seq_name}->{$chr}->fill($hole_length);
        }

        for my $set ( $info_of{$seq_name}->{$chr}->sets ) {
            my $piece_coverage = $set->size / $length;
            my $location = sprintf "%s:%s-%s", $chr, $set->min, $set->max;
            if ( $piece_coverage >= $coverage ) {
                push @locations, $location;
                $count_of{$real_name}++;
            }
        }
    }
}

for my $seq_name ( keys %count_of ) {
    if ( $count_of{$seq_name} != 1 ) {
        print "    $seq_name got $count_of{$seq_name}\n";
    }
}

$stopwatch->block_message("Retrieve sequences from genome");
@locations = uniq(@locations);
for my $location (@locations) {
    my $seq_genome = get_seq_faidx( $genome, $location );
    path($output)->append(">$location\n");
    path($output)->append("$seq_genome\n");
}

if ($debug) {
    $stopwatch->block_message("Dump mem.yml");
    DumpFile( 'mem.yml',
        { mem_result => $mem_result, info => \%info_of, locations => \@locations } );
}

$stopwatch->end_message;

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub run_sparsemem {
    my $file   = shift;
    my $genome = shift;
    my $length = shift || 20;

    my $cmd = sprintf "sparsemem -maxmatch -F -l %d -b -n -k 3 -threads 3 %s %s", $length, $genome,
        $file;
    my $result = `$cmd`;

    return $result;
}

__END__
