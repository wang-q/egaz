#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Basename;
use IO::Zlib;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_file;    # blocked fasta file
my $out_dir;         # Specify output dir here

my $gzip;                     # open .gz

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'i|in_file=s'  => \$in_file,
    'o|out_dir=s'  => \$out_dir,
    'gzip'         => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;

# make output dir
unless ($out_dir) {
    $out_dir = File::Spec->rel2abs($in_file) . "_vcf";
}
if ( !-e $out_dir ) {
    mkdir $out_dir, 0777;
}
else {
    print "We are going to output blocked fasta.\n";
    die "$out_dir exists, you should remove it first to avoid errors.\n";
}

my $in_fh;
if ( !$gzip ) {
    open $in_fh, '<', $in_file;
}
else {
    $in_fh = IO::Zlib->new( $in_file, "rb" );
}

{
    my $content = '';
    while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';
            die "headers not equal to seqs\n" if @lines % 2;
            die "Two few lines in block\n" if @lines < 4;

            my ( @headers, %seq_of, @simple_names );
            while (@lines) {
                # header
                my $header = shift @lines;
                $header =~ s/^\>//;
                chomp $header;
                push @headers, $header;
                
                # seq
                my $seq = shift @lines;
                chomp $seq;
                $seq = uc $seq;
                $seq_of{$header} = $seq;
                
                # simple name
                my $info_ref = decode_header($header);
                push @simple_names, $info_ref->{name};
            }
            
            # build info required by write_fasta()
            my $out_filename = $headers[0];
            $out_filename =~ s/\|.+//; # remove addtional fields
            $out_filename =~ s/[\(\)\:]+/./g;
            $out_filename .= '.fas';
            my $out_file = File::Spec->catfile($out_dir, $out_filename);
            write_fasta($out_file, \%seq_of, \@headers, \@simple_names);
        }
        else {
            $content .= $line;
        }
    }
}

if ( !$gzip ) {
    close $in_fh;
}
else {
    $in_fh->close;
}

$stopwatch->block_message( "Process finished.", "duration" );

exit;

__END__

=head1 NAME

    split_fas.pl - split blocked fasta file to seperate per-alignment files.

=head1 SYNOPSIS

    perl split_fas.pl --in I.net.axt.fas 

      Options:
        --help              brief help message
        --man               full documentation
        -i, --in_file       blocked fasta file's location
        -o, --out_dir       output location
        --gzip              input file is gzipped

=cut

perl split_fas.pl --in ~/data/alignment/yeast_combine/S288CvsIII_mft/chrI.synNet.maf.gz.fas 