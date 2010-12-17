#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use File::Find::Rule;
use Path::Class;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_dir;
my $out_dir;

my $chunk_size = 10_010_000;
my $overlap    = 10_000;

my $wrap_length = 50;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    'in_dir=s'      => \$in_dir,
    'out_dir=s'     => \$out_dir,
    'chunk_size=i'  => \$chunk_size,
    'overlap=i'     => \$overlap,
    'wrap_length=i' => \$wrap_length,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# now run!
#----------------------------------------------------------#

my @files = File::Find::Rule->file->name('*.fa')->in($in_dir);
die "You should provide a dir or file.\n" unless @files;

# make output dir
unless ( -e $out_dir ) {
    mkdir $out_dir, 0777
        or die "Cannot create \"$out_dir\" directory: $!";
}

my %length_of;
for my $file (@files) {
    print "Read out file: $file\n";
    my ( $seq_of, $seq_names ) = read_fasta($file);
    for my $name ( @{$seq_names} ) {
        print "Fasta header: $name\n";
        if ( exists $length_of{$name} ) {
            print "Seq $name has been processed, skip it.\n";
            next;
        }

        my $new_name = $name;
        if ( $name =~ /[^\w]/ ) {
            $new_name =~ s/[^\w].*//;
            print "$name is too complicated, simplify it to $new_name\n";
        }

        write_seq( $out_dir, $new_name, $seq_of->{$name}, $wrap_length );

        my $size = length $seq_of->{$name};
        if ( $size > $chunk_size + $overlap ) {

            # break it up
            my $intervalsRef
                = overlappingIntervals( 0, $size, $chunk_size, $overlap );
            for my $i ( @{$intervalsRef} ) {
                my ( $start, $end ) = @{$i};
                write_empty_file( $out_dir, $new_name, $start, $end );
            }

        }
        else {
            write_empty_file( $out_dir, $new_name, 1, $size );
        }
        $length_of{$new_name} = $size;
    }
}

{
    my $file = file( $out_dir, "length.csv" );
    my $fh = $file->openw;
    for my $key ( sort keys %length_of ) {
        print {$fh} "$key,$length_of{$key}\n";
    }
}

sub overlappingIntervals {

   # Return a list of overlapping intervals (1-based starts if $opt_oneBased).
   # Starts increment by $chunk; each end overlaps the next start by $lap.
    my ( $start, $end, $chunk, $lap ) = @_;
    my @intervals;
    for ( my $iStart = $start; $iStart < $end; $iStart += $chunk ) {
        my $iEnd = $iStart + $chunk + $lap;
        $iEnd = $end if ( $iEnd > $end );
        push @intervals, [ ( $iStart + 1 ), $iEnd ];
    }
    return \@intervals;
}

sub write_seq {
    my $dir         = shift;
    my $name        = shift;
    my $seq         = shift;
    my $wrap_length = shift;

    my $seq_length = length $seq;

    my $file = file( $dir, "$name.fa" );
    my $fh = $file->openw;
    print {$fh} ">$name\n";

    if ($wrap_length) {
        for ( my $pos = 0; $pos < $seq_length; $pos += $wrap_length ) {
            print {$fh} substr( $seq, $pos, $wrap_length ), "\n";
        }
    }
    else {
        print {$fh} $seq, "\n";
    }

}

sub write_empty_file {
    my $dir   = shift;
    my $name  = shift;
    my $start = shift;
    my $end   = shift;

    my $file = file( $dir, "$name.fa[$start,$end]" );
    my $fh = $file->openw;

    return;
}

exit;

=head1 NAME

    part_seq.pl - Partitions a set of sequences by size

=head1 SYNOPSIS

    perl part_seq.pl -in t\S288C -out t\S288C_parted -chunk 500000
    perl part_seq.pl -in t\RM11 -out t\parted -chunk 500000

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Blastz will take the first sequence in target fasta file and all sequences in
query fasta file.

So if there are mutiple query files, the program uses the largest one. And all
target files in dir_target will be processed.

So, if there are combined fasta files and multi fasta files coexisting in the
target directory, just delete the axt file matched with the combined fasta
filename.

=cut
