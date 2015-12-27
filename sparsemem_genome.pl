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
        --length            INT     Match length of sparsemem, default is [50]
        --coverage      -c  FLOAT   default is [0.9]
        --hole              INT     Fill holes (snp/indel) less than this length, default is [10]
        --debug                     Write mem.yml for debugging

=head1 REQUIREMENTS

    brew tap homebrew/science
    brew tap wang-q/tap
    brew install sparsemem samtools

=cut

my $output;

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \my $file,
    'genome|g=s' => \my $genome,
    'coverage|c=f' => \( my $coverage     = 0.9 ),
    'hole=i'       => \( my $hole_length  = 10 ),
    'length=i'     => \( my $match_length = 50 ),
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
    close($fh_pipe);

    return $seq;
}

sub decode_header {
    my $header = shift;

    # S288C.chrI(+):27070-29557|species=S288C
    my $head_qr = qr{
                ([\w_]+)?           # name
                [\.]?               # spacer
                ((?:chr)?[\w-]+)    # chr name
                (?:\((.+)\))?       # strand
                [\:]                # spacer
                (\d+)               # chr start
                [\_\-]              # spacer
                (\d+)               # chr end
            }xi;

    tie my %info, "Tie::IxHash";

    $header =~ $head_qr;
    my $name     = $1;
    my $chr_name = $2;

    if ( defined $name ) {
        %info = (
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        );
        if ( !defined $info{chr_strand} ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '1' ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '-1' ) {
            $info{chr_strand} = '-';
        }
    }
    elsif ( defined $chr_name ) {
        $name = $header;
        %info = (
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        );
        if ( !defined $info{chr_strand} ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '1' ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '-1' ) {
            $info{chr_strand} = '-';
        }
    }
    else {
        $name = $header;
        %info = (
            chr_name   => 'chrUn',
            chr_strand => '+',
            chr_start  => undef,
            chr_end    => undef,
        );
    }
    $info{name} = $name;

    # additional keys
    if ( $header =~ /\|(.+)/ ) {
        my @parts = grep {defined} split /;/, $1;
        for my $part (@parts) {
            my ( $key, $value ) = split /=/, $part;
            if ( defined $key and defined $value ) {
                $info{$key} = $value;
            }
        }
    }

    return \%info;
}

__END__
