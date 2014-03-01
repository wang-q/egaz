#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Class;
use List::Util qw(reduce);
use List::MoreUtils qw(any all uniq zip);
use String::Compare;
use Set::Scalar;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $size_file;

my $file;

my $multi_key;    # file has multiple keys

my $outfile;

my $remove_chr;    # remove "chr" in "chr1"

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    's|size=s'    => \$size_file,
    'f|file=s'    => \$file,
    'mk'          => \$multi_key,
    'o|outfile=s' => \$outfile,
    'r|remove'    => \$remove_chr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

if ( !$outfile ) {
    $outfile = "$file.csv";
}

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Compare runlist...");

#----------------------------#
# Loading
#----------------------------#
my $length_of = read_sizes( $size_file, $remove_chr );

print "Loading $file...\n";
my $s_of = {};
my @keys;
if ($multi_key) {
    my $yml = LoadFile($file);
    @keys = sort keys %{$yml};

    for my $key (@keys) {
        $s_of->{$key} = runlist2set( $yml->{$key}, $remove_chr );
    }
}
else {
    @keys = ("__single");
    $s_of->{__single} = runlist2set( LoadFile($file), $remove_chr );
}

#----------------------------#
# Calcing
#----------------------------#
print "\nCalc and write to $outfile\n";

if ($multi_key) {
    my @lines;

    for my $key (@keys) {
        my @key_lines = csv_lines( $s_of->{$key}, $length_of );
        $_ = "$key,$_" for @key_lines;
        push @lines, @key_lines;
    }

    unshift @lines, "key,name,length,size,coverage\n";
    open my $fh, '>', $outfile;
    print {$fh} $_ for @lines;
    close $fh;
}
else {
    my @lines = csv_lines( $s_of->{__single}, $length_of );

    unshift @lines, "name,length,size,coverage\n";
    open my $fh, '>', $outfile;
    print {$fh} $_ for @lines;
    close $fh;
}

$stopwatch->end_message;

sub read_sizes {
    my $file       = shift;
    my $remove_chr = shift;

    my $fh = file($file)->openr;
    my %length_of;
    while (<$fh>) {
        chomp;
        my ( $chr, $value ) = split /\t/;
        $chr =~ s/chr0?// if $remove_chr;
        $length_of{$chr} = $value;
    }

    return \%length_of;
}

sub runlist2set {
    my $runlist_of = shift;
    my $remove_chr = shift;

    my $set_of = {};

    for my $chr ( sort keys %{$runlist_of} ) {
        my $new_chr = $chr;
        $new_chr =~ s/chr0?// if $remove_chr;
        my $set = AlignDB::IntSpan->new( $runlist_of->{$chr} );
        printf "\tkey:\t%s\tlength:\t%s\n", $new_chr, $set->size;
        $set_of->{$new_chr} = $set;
    }

    return $set_of;
}

sub csv_lines {
    my $set_of    = shift;
    my $length_of = shift;

    my @lines;

    my ( $all_length, $all_size, $all_coverage );
    for my $chr ( sort keys %{$set_of} ) {
        my $length   = $length_of->{$chr};
        my $size     = $set_of->{$chr}->size;
        my $coverage = $size / $length;

        $all_length += $length;
        $all_size   += $size;

        push @lines, "$chr,$length,$size,$coverage\n";
    }

    $all_coverage = $all_size / $all_length;
    push @lines, "all,$all_length,$all_size,$all_coverage\n";

    return @lines;
}

__END__

=head1 NAME

    stat_runlist.pl - coverage on chromosomes for runlists

=head1 SYNOPSIS

    perl stat_runlist.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --op                operations: intersect, union, diff
        --file1
        --file2
        --outfile

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

