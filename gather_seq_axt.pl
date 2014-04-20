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

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $dir;
my $length = 1000;
my $gzip;    # open .gz

my $avoid_regex;    # '(random|hap|Un|M)';

my $output;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'd|dir=s'    => \$dir,
    'l|length=s' => \$length,
    'g|gzip'     => \$gzip,
    'r|regex=s'  => \$avoid_regex,
    'o|output=s' => \$output,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Merge chr runlist...");

if ( !$output ) {
    $output = basename($dir);
    ($output) = grep {defined} split /\./, $output;
}

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt')->in($dir);
    printf "\n----Total .axt Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt.gz')->in($dir);
    printf "\n----Total .axt.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
open my $fasta_fh, ">", "$output.axt.fasta";
for my $file (@files) {
    $stopwatch->block_message("Analysis [$file]");
    my $data = parse_axt( $file, 1 );
    @{$data} = grep { $_->[2] >= $length } @{$data};

    if ($avoid_regex) {
        @{$data} = grep { $_->[0]{chr_name} !~ /$avoid_regex/i } @{$data};
        @{$data} = grep { $_->[1]{chr_name} !~ /$avoid_regex/i } @{$data};
    }

    for my $info_ref ( @{$data} ) {
        for my $i ( 0, 1 ) {
            my $chr_name   = $info_ref->[$i]{chr_name};
            my $chr_strand = $info_ref->[$i]{chr_strand};
            my $chr_start  = $info_ref->[$i]{chr_start};
            my $chr_end    = $info_ref->[$i]{chr_end};

            # piece name
            my $head = "$chr_name:$chr_start-$chr_end";

            my $seq = $info_ref->[$i]{seq};
            $seq =~ tr/-//d;
            $seq = uc $seq;
            if ($chr_strand ne "+") {
                $seq = revcom($seq);
            }
            print {$fasta_fh} ">$head\n";
            print {$fasta_fh} "$seq\n";
        }
    }
}
close $fasta_fh;

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

sub parse_axt {
    my $file = shift;
    my $gzip = shift;

    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $file;
    }
    else {
        $in_fh = IO::Zlib->new( $file, "rb" );
    }

    my @data;
    while (1) {
        my $summary_line = <$in_fh>;
        last unless $summary_line;
        next if $summary_line =~ /^#/;

        chomp $summary_line;
        chomp( my $first_line  = <$in_fh> );
        chomp( my $second_line = <$in_fh> );
        my $dummy = <$in_fh>;    # blank line

        my ($align_serial, $first_chr,    $first_start,
            $first_end,    $second_chr,   $second_start,
            $second_end,   $query_strand, $align_score,
        ) = split /\s+/, $summary_line;

        my $info_refs = [
            {   chr_name   => $first_chr,
                chr_start  => $first_start,
                chr_end    => $first_end,
                chr_strand => '+',
                seq        => $first_line,
            },
            {   chr_name   => $second_chr,
                chr_start  => $second_start,
                chr_end    => $second_end,
                chr_strand => $query_strand,
                seq        => $second_line,
            },
            length $first_line,
        ];

        push @data, $info_refs;
    }

    if ( !$gzip ) {
        close $in_fh;
    }
    else {
        $in_fh->close;
    }

    return \@data;
}

sub revcom {
    my $seq = shift;

    $seq =~ tr/ACGTMRWSYKVHDBNacgtmrwsykvhdbn-/TGCAKYWSRMBDHVNtgcakyswrmbdhvn-/;
    my $seq_rc = reverse $seq;

    return $seq_rc;
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

    perl gather_info_axt.pl -d d:\data\alignment\self_alignment\S288Cvsselfalign\

=cut


