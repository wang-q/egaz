#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use IO::Zlib;
use Path::Tiny;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(read_sizes);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

fas2vcf.pl - list variations in blocked fasta file

=head1 SYNOPSIS

    perl fas2vcf.pl --in <.fas> -s <chr.sizes> [options]
      Options:
        --help              brief help message
        --man               full documentation
        -i, --in_file       blocked fasta file's location
        -o, --out_file      output location
        -s, --size_file     chr.size

    perl fas2vcf.pl --in I.net.axt.fas -s chr.sizes

=cut

my $jvk = '~/share/jvarkit/biostar94573.jar';

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'i|in_file=s'   => \my $in_file,
    'o|out_file=s'  => \my $out_file,
    's|size_file=s' => \my $size_file,
) or HelpMessage(1);

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;

if ( !defined $out_file ) {
    $out_file = path($in_file)->absolute->stringify . ".vcf";
}

my $temp_dir = Path::Tiny->tempdir( TEMPLATE => "fas_XXXXXXXX" );
my $temp_list = Path::Tiny->tempfile( TEMPLATE => "fas_XXXXXXXX" );
my $length_of = read_sizes($size_file);

{    #split fas
    my $cmd = "fasops split $in_file --rm -o $temp_dir";
    exec_cmd($cmd);
}

{
    my @files = $temp_dir->children(qr/\.fas$/);

    for my $f (@files) {
        my ( $name, $chr_name, $chr_strand, $chr_pos ) = split /\./, $f->basename(".fas");
        my ( $chr_start, $chr_end ) = split /\-/, $chr_pos;
        my $chr_length = $length_of->{$chr_name};

        my $cmd = "java -jar $jvk $f";
        my @lines = split /\n/, `$cmd`;

        for my $l (@lines) {
            if ( $l =~ /^\#\#contig\=\<ID\=/ ) {
                $l = "##contig=<ID=$chr_name,length=@{[$length_of->{$chr_name}]}>";
            }
            if ( $l =~ /^chrUn\t/ ) {
                my @fields = split /\t/, $l;
                $fields[0] = $chr_name;

                # vcf position is 1-based
                $fields[1] = $chr_start + $fields[1] - 1;
                $l = join "\t", @fields;
            }
        }

        my $vf = $f . '.vcf';
        path($vf)->spew( map { $_ . "\n" } @lines );
        $temp_list->append( $vf . "\n" );
    }
}

{    # concat vcf
    my $cmd = "vcf-concat -f $temp_list > $out_file";
    exec_cmd($cmd);
}

$stopwatch->block_message( "Process finished.", "duration" );

exit;

sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

__END__
