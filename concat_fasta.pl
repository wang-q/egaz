#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use MCE;

use File::Find::Rule;
use IO::Zlib;
use Path::Tiny;
use List::MoreUtils qw(uniq);

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

concat_fasta.pl - concatenate blocked fasta files

=head1 SYNOPSIS

    perl concat_fasta.pl [options]
      Options:
        --help          -?          brief help message
        --in_dir        -i  STR     fasta files' location
        --out_file      -o  STR     output file
        --sampling      -s          random sampling
        --total_length  -l  INT     when exceed this, quit. Default is [1_000_000]
        --relaxed       -rp         output relaxed phylip instead of fasta
        --wrap          -w  INT     long line wrap length
        --gzip                      input files are gzipped. .fas.gz

    perl concat_fasta.pl --in_dir S288CvsRM11 -o S288CvsRM11.fasta

=cut

GetOptions(
    'help|?'       => sub { HelpMessage(0) },
    'in_dir|i=s'   => \my $in_dir,
    'out_file|o=s' => \my $out_file,
    'sampling|s'   => \my $sampling,
    'total|l=i' => \( my $total_length = 1_000_000 ),
    'relaxed|rp' => \my $relaxed_phylip,
    'wrap|w=i'   => \( my $wrap = 0 ),
    'gzip'       => \my $gzip,
) or HelpMessage(1);

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
unless ($out_file) {
    $out_file = path($in_dir)->stringify . ( $relaxed_phylip ? ".concat.phy" : ".concat.fas" );
}

my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($in_dir);
    printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name( '*.fa.gz', '*.fas.gz', '*.fasta.gz' )->in($in_dir);
    printf "\n----Total .fas.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

exit unless scalar @files;

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
# process each .fasta files
my $stopwatch = AlignDB::Stopwatch->new;

my $all_names = gather_names( \@files, $gzip );
print Dump $all_names;

my $seq_of_ary = [];
for my $file (@files) {
    print "Process $file\n";
    my $count = gather_seq( $file, $all_names, $seq_of_ary, $gzip );
    print "Add $count sequence blocks\n";
}

my $all_seq_of = {};
if ($sampling) {
    while (1) {
        my $rand_idx = int rand $#{$seq_of_ary};
        $all_seq_of->{$_} .= $seq_of_ary->[$rand_idx]{$_} for @{$all_names};
        last if length $all_seq_of->{ $all_names->[0] } > $total_length;
    }
}
else {
    for my $idx ( 0 .. $#{$seq_of_ary} ) {
        $all_seq_of->{$_} .= $seq_of_ary->[$idx]{$_} for @{$all_names};
    }
}
my $seq_length = length $all_seq_of->{ $all_names->[0] };

open my $out_fh, '>', $out_file;
if ($relaxed_phylip) {
    print {$out_fh} scalar @{$all_names}, " $seq_length\n";
    for my $name ( @{$all_names} ) {
        print {$out_fh} "$name ";
        print {$out_fh} $all_seq_of->{$name}, "\n";
    }
}
else {
    for my $name ( @{$all_names} ) {
        print {$out_fh} ">$name\n";
        if ($wrap) {
            for ( my $pos = 0; $pos < $seq_length; $pos += $wrap ) {
                print {$out_fh} substr( $all_seq_of->{$name}, $pos, $wrap ), "\n";
            }
        }
        else {
            print {$out_fh} $all_seq_of->{$name}, "\n";
        }
    }
}
close $out_fh;

$stopwatch->block_message( "All files have been processed.", "duration" );
exit;

sub gather_names {
    my $files = shift;
    my $gzip  = shift;

    my @names;
    for my $file ( @{$files} ) {
        my $in_fh;
        if ( !$gzip ) {
            open $in_fh, '<', $file;
        }
        else {
            $in_fh = IO::Zlib->new( $file, "rb" );
        }
        while (<$in_fh>) {
            my $index = index $_, '>';
            if ( $index == 0 ) {
                chomp;
                push @names, $_;
            }
        }
        close $in_fh;
    }

    @names = map { s/\>//; ( split /\s+/ )[0] } @names;
    @names = map { ( split /\./ )[0] } @names;
    @names = uniq(@names);

    return \@names;
}

sub gather_seq {
    my $infile     = shift;
    my $all_names  = shift;
    my $ary_seq_of = shift;
    my $gzip       = shift;

    my $count = 0;

    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $infile;
    }
    else {
        $in_fh = IO::Zlib->new( $infile, "rb" );
    }

    my $content = '';
    while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';
            die "headers not equal to seqs\n" if @lines % 2;
            die "Two few lines in block\n" if @lines < 4;

            my ( $seq_of, $seq_names ) = ( {}, [] );
            while (@lines) {
                my $name = shift @lines;
                chomp $name;
                $name =~ s/^\>//;
                $name = ( split /\s+/, $name )[0];
                $name = ( split /\./,  $name )[0];
                push @{$seq_names}, $name;

                my $seq = shift @lines;
                chomp $seq;
                $seq_of->{$name} = uc $seq;
            }

            my $align_length = length $seq_of->{ $seq_names->[0] };
            for my $name ( @{$seq_names} ) {
                if ( ( length $seq_of->{$name} ) != $align_length ) {
                    die "Sequences should have the same length!\n";
                }
            }

            # concat seqs
            # fill unpresented names with ------
            for my $name ( @{$all_names} ) {
                my $flag = grep { $_ eq $name } @{$seq_names};
                if ( !$flag ) {
                    $seq_of->{$name} = '-' x $align_length;
                }
            }

            # collect seq_of
            push @{$seq_of_ary}, $seq_of;

            $count++;
        }
        else {
            $content .= $line;
        }
    }
    close $in_fh;

    return $count;
}

__END__
