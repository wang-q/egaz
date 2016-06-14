#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use IO::Zlib;
use Path::Tiny;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

fas2vcf.pl - list variations in blocked fasta file

=head1 SYNOPSIS

    perl fas2vcf.pl -i <.fas> -s <chr.sizes> [options]
      Options:
        --help          -?          brief help message
        --input         -i  STR     input blocked fasta file
        --output        -o  STR     output vcf file
        --size          -s  STR     chr.size

    perl fas2vcf.pl -i I.net.axt.fas -s chr.sizes

=cut

my $jvk = '~/share/jvarkit/biostar94573.jar';

GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'input|i=s'  => \my $in_file,
    'output|o=s' => \my $out_file,
    'size|s=s'   => \my $size_file,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;

if ( !defined $out_file ) {
    $out_file = path($in_file)->absolute->stringify . ".vcf";
}

my $temp_dir = Path::Tiny->tempdir( TEMPLATE => "fas_XXXXXXXX" );
my $temp_list = Path::Tiny->tempfile( TEMPLATE => "fas_XXXXXXXX" );
my $length_of = App::RL::Common::read_sizes($size_file);

{
    $stopwatch->block_message("Split fas");
    my $cmd = "fasops split $in_file --rm -o $temp_dir";
    exec_cmd($cmd);
}

{
    $stopwatch->block_message("Fas2vcf");
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

{
    $stopwatch->block_message("concat vaf");
    my $cmd = "vcf-concat -f $temp_list > $out_file";
    exec_cmd($cmd);
}

$stopwatch->block_message( "Process finished.", "duration" );

exit;

__END__
