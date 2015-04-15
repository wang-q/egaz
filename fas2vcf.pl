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
use Path::Tiny;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_file;     # blocked fasta file
my $out_file;    # Specify output dir here

my $size_file;   # chr.sizes
my $jvk = '~/share/jvarkit/biostar94573.jar';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    'i|in_file=s'   => \$in_file,
    'o|out_file=s'  => \$out_file,
    's|size_file=s' => \$size_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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
    my @files = $temp_dir->children( qr/\.fas$/ );

    for my $f (@files) {
        my ( $name, $chr_name, $chr_strand, $chr_pos ) = split /\./,
            $f->basename(".fas");
        my ( $chr_start, $chr_end ) = split /\-/, $chr_pos;
        my $chr_length = $length_of->{$chr_name};

        my $cmd = "java -jar $jvk $f";
        my @lines = split /\n/, `$cmd`;

        for my $l (@lines) {
            if ( $l =~ /^\#\#contig\=\<ID\=/ ) {
                $l
                    = "##contig=<ID=$chr_name,length=@{[$length_of->{$chr_name}]}>";
            }
            if ( $l =~ /^chrUn\t/ ) {
                my @fields = split /\t/, $l;
                $fields[0] = $chr_name;
                $fields[1]
                    = $chr_start + $fields[1] - 1;    # vcf position is 1-based
                $l = join "\t", @fields;
            }
        }

        my $vf = $f . '.vcf';
        path($vf)->spew(map {$_ . "\n"} @lines);
        $temp_list->append($vf . "\n");
    }
}

{ # concat vcf
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

=head1 NAME

    fas2vcf.pl - list variations in blocked fasta file

=head1 SYNOPSIS

    perl fas2vcf.pl --in I.net.axt.fas -s chr.sizes

      Options:
        --help              brief help message
        --man               full documentation
        -i, --in_file       blocked fasta file's location
        -o, --out_file      output location
        -s, --size_file     chr size

=cut

perl fas2vcf.pl -i ~/data/alignment/yeast_combine/S288CvsIII_mft/chrI.synNet.maf.gz.fas \
    -s /Users/wangq/data/alignment/yeast_combine/S288C/chr.sizes \
    
