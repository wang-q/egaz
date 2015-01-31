#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Find::Rule;
use File::Basename;
use IO::Zlib;

use MCE;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_dir = '.';    # Specify location here
my $out_dir;         # Specify output dir here
my $length = 1000;   # Set the threshold of alignment length

my $tname = 'target';
my $qname = 'query';

# run in parallel mode
my $parallel = 1;

my $gzip;            # open .gz

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'i|in_dir=s'  => \$in_dir,
    'o|out_dir=s' => \$out_dir,
    'l|length=i'  => \$length,
    't|tname=s'   => \$tname,
    'q|qname=s'   => \$qname,
    'parallel=i'  => \$parallel,
    'gzip'        => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
unless ($out_dir) {
    $out_dir = File::Spec->rel2abs($in_dir) . "_fasta";
}
if ( !-e $out_dir ) {
    mkdir $out_dir, 0777;
}

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt')->in($in_dir);
    printf "\n----Total .axt Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name('*.axt.gz')->in($in_dir);
    printf "\n----Total .axt.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    # open files
    my $infile = $chunk_ref->[0];
    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $infile;
    }
    else {
        $in_fh = IO::Zlib->new( $infile, "rb" );
    }

    my $outfile = basename($infile);
    $outfile = $out_dir . "/$outfile" . ".fas";
    open my $out_fh, ">", $outfile;

    # read and write
    my $data = parse_axt( $infile, $gzip );
    @{$data} = grep { $_->[2] >= $length } @{$data};

    for my $info_ref ( @{$data} ) {
        for my $i ( 0, 1 ) {
            $info_ref->[$i]{name} = $i == 0 ? $tname : $qname;
            printf {$out_fh} ">%s\n", encode_header( $info_ref->[$i] );
            printf {$out_fh} "%s\n",  $info_ref->[$i]{seq};
        }
        print {$out_fh} "\n";
    }

    # close file handlers
    if ( !$gzip ) {
        close $in_fh;
    }
    else {
        $in_fh->close;
    }
    close $out_fh;

    print ".fas file generated.\n\n";
};

# process each files
my $stopwatch = AlignDB::Stopwatch->new;

my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
$mce->foreach( [ sort @files ], $worker );

$stopwatch->block_message( "All files have been processed.", "duration" );

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

__END__


=head1 NAME

    axt2fasta.pl - convert axt to fasta

=head1 SYNOPSIS
    perl axt2fasta.pl --in G:/S288CvsRM11

    axt2fasta.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --in                axt files' location
        --out               output location
        --length            length threshold
        --parallel          run in parallel mode

=cut

perl axt2fasta.pl --in ~/data/alignment/self/yeast_new/S288cvsselfalign \
    -tname S288c -qname S288c

# Merge all fas files
find ~/data/alignment/self/yeast_new/S288cvsselfalign -name "*.fas" \
    | sort \
    | xargs perl -nl -e '/^>/ and print and next; /^$/ and next; $_ = uc $_; $_ =~ tr/-//d; print'
