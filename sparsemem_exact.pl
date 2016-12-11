#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use Tie::IxHash;
use List::MoreUtils qw(uniq);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use App::RL::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

sparsemem_exact.pl - Get exact genome locations of sequence pieces in a fasta file

=head1 SYNOPSIS

    perl sparsemem_exact.pl -f <blast result file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     query fasta file
        --genome        -g  STR     reference genome file
        --length        -l  INT     Match length of sparsemem, default is [100]
        --discard       -d  INT     Discard blocks with copy number larger than this,
                                    default is [0]
        --output        -o  STR     output
        --debug                     Write mem.yml for debugging

=head1 REQUIREMENTS

    brew tap homebrew/science
    brew tap wang-q/tap
    brew install sparsemem samtools

=cut

GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $file,
    'genome|g=s' => \my $genome,
    'length|l=i'  => \( my $match_length = 100 ),
    'discard|d=i' => \( my $discard      = 0 ),
    'output|o=s'  => \( my $output ),
    'debug'       => \( my $debug ),
) or Getopt::Long::HelpMessage(1);

if ( !$output ) {
    $output = path($file)->basename( ".fasta", ".fas", ".fa" );
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.replace.txt";
}
path($output)->remove;

if ( !defined $file ) {
    die "Need --file\n";
}
elsif ( !path($file)->is_file ) {
    die "--file [$file] don't exist\n";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Get locations...");

#----------------------------------------------------------#
# load blast reports
#----------------------------------------------------------#
$stopwatch->block_message("Generate sparsemem reports");

my $mem_result = run_sparsemem( $file, $genome, $match_length );

$stopwatch->block_message("Parse sparsemem");
tie my %match_of, "Tie::IxHash";
my $fh_mem = $mem_result->openr;
my $cur_name;
while ( my $line = <$fh_mem> ) {
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

        if ( !exists $match_of{$cur_name} ) {
            $match_of{$cur_name} = [];
        }
        push @{ $match_of{$cur_name} }, \@fields;
    }
    else {    # Blank line, do nothing
    }
}
close $fh_mem;

$stopwatch->block_message("Get exact matches");
my $size_of = get_size_faops($file);
tie my %locations, "Tie::IxHash";    #  genome location - simple location
for my $seq_name ( keys %match_of ) {
    my $ori_name   = $seq_name;
    my $chr_strand = "+";
    if ( $ori_name =~ / Reverse/ ) {
        $ori_name =~ s/ Reverse//;
        $chr_strand = "-";
    }
    my $species_name = App::RL::Common::decode_header($ori_name)->{name};

    for my $f ( @{ $match_of{$seq_name} } ) {
        next unless $f->[2] == 1;
        next unless $f->[3] == $size_of->{$ori_name};

        my $header = App::RL::Common::encode_header(
            {   name   => $species_name,
                chr    => $f->[0],
                start  => $f->[1],
                end    => $f->[1] + $f->[3] - 1,
                strand => $chr_strand,
            }
        );
        if ( !exists $locations{$ori_name} ) {
            $locations{$ori_name} = {};
        }

        $locations{$ori_name}->{$header} = 1;
    }
}

$stopwatch->block_message("Write replace tsv file");
for my $ori_name ( keys %locations ) {
    my @matches = keys %{ $locations{$ori_name} };
    if ( $discard and @matches > $discard ) {
        printf "    %s\tgot %d matches, discard it\n", $ori_name, scalar(@matches);
        path($output)->append($ori_name);
        path($output)->append("\n");
    }
    elsif ( @matches > 1 ) {
        printf "    %s\tgot %d matches\n", $ori_name, scalar(@matches);
        path($output)->append( join( "\t", $ori_name, @matches ) );
        path($output)->append("\n");
    }
    elsif ( $ori_name ne $matches[0] ) {
        printf "    %s\tgot wrong position, %s\n", $ori_name, $matches[0];
        path($output)->append( join( "\t", $ori_name, @matches ) );
        path($output)->append("\n");
    }
}

if ($debug) {
    $stopwatch->block_message("Dump mem_exact.yml");
    YAML::Syck::DumpFile(
        'mem_exact.yml',
        {   mem_result => $mem_result,
            info       => \%match_of,
            locations  => \%locations
        }
    );
}

$stopwatch->end_message;

exit;

sub run_sparsemem {
    my $file   = shift;
    my $genome = shift;
    my $length = shift || 20;

    my $result = Path::Tiny->tempfile;

    my $cmd = sprintf "sparsemem -maxmatch -F -l %d -b -n -k 3 -threads 3 %s %s > %s",
        $length,
        $genome,
        $file,
        $result->stringify;
    system $cmd;

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
    close $fh_pipe;

    return $seq;
}

sub get_size_faops {
    my $file = shift;

    my $cmd = sprintf "faops size %s", $file;
    my @lines = grep {defined} split /\n/, `$cmd`;

    tie my %length_of, "Tie::IxHash";
    for (@lines) {
        my ( $key, $value ) = split /\t/;
        $length_of{$key} = $value;
    }

    return \%length_of;
}

__END__
