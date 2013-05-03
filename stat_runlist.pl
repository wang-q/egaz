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

my $outfile;

my $remove_chr;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    's|size=s'    => \$size_file,
    'f|file=s'    => \$file,
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

my $length_of = read_sizes( $size_file, $remove_chr );

print "Loading...\n";
my $set_of = {};
{
    print "Filename: $file\n";
    my $runlist_of = LoadFile($file);

    for my $key ( sort keys %{$runlist_of} ) {
        my $set = AlignDB::IntSpan->new( $runlist_of->{$key} );
        $key =~ s/chr0?// if $remove_chr;
        printf "key:\t%s\tlength:\t%s\n", $key, $set->size;
        $set_of->{$key} = $set;
    }
}
print "\n";

open my $fh, '>', $outfile;
print {$fh} "name,length,size,coverage\n";
for my $key ( sort keys %{$set_of} ) {
    print "For $key\n";
    if ( exists $length_of->{$key} ) {
        my $length   = $length_of->{$key};
        my $size     = $set_of->{$key}->size;
        my $coverage = $size / $length;
        print "$key,$length,$size,$coverage\n";
        print {$fh} "$key,$length,$size,$coverage\n";
    }
}
close $fh;

$stopwatch->end_message;

sub read_sizes {
    my $file       = shift;
    my $remove_chr = shift;

    my $fh = file($file)->openr;
    my %length_of;
    while (<$fh>) {
        chomp;
        my ( $key, $value ) = split /\t/;
        $key =~ s/chr0?// if $remove_chr;
        $length_of{$key} = $value;
    }

    return \%length_of;
}

sub read_bed {
    my $file       = shift;
    my $remove_chr = shift;
    my @data;

    open my $fh, '<', $file;
    while ( my $string = <$fh> ) {
        next unless defined $string;
        chomp $string;
        my ( $chr, $start, $end )
            = ( split /\t/, $string )[ 0, 1, 2 ];
        next unless $chr =~ /^\w+$/;
        $chr =~ s/chr0?//i if $remove_chr;
        next unless $start =~ /^\d+$/;
        next unless $end =~ /^\d+$/;
        if ( $start > $end ) {
            ( $start, $end ) = ( $end, $start );
        }
        next if $end - $start < 10;
        my $set = AlignDB::IntSpan->new("$start-$end");
        push @data, { chr => $chr, set => $set, };
    }
    close $fh;

    return \@data;
}

__END__

=head1 NAME

    compare_runlist.pl - compare 2 chromosome runlists

=head1 SYNOPSIS

    perl compare_runlist.pl [options]
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

