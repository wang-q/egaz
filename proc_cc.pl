#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use File::Spec;
use Path::Class;
use Graph;
use List::Util qw(sum0);
use List::MoreUtils qw(zip);

use Bio::SeqIO;

use FindBin;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $cc_file;
my $size_file;

my $output;

my $low_cut = 4;    # low copies their own stats

my $write_seq;

my $verbose;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'f|file=s'    => \$cc_file,
    's|size=s'    => \$size_file,
    'o|output=s'  => \$output,
    'l|low=s'     => \$low_cut,
    'w|write_seq' => \$write_seq,
    'v|verbose'   => \$verbose,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$cc_file]");

if ( !$output ) {
    $output = basename($cc_file);

    ($output) = grep {defined} split /\./, $output;
    $output = "$output.cc";
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

#----------------------------#
# load cc
#----------------------------#
my $yml = LoadFile($cc_file);

my @cc         = @{ $yml->{cc} };
my $count_of   = $yml->{count};
my $runlist_of = $yml->{runlist};

#----------------------------#
# coverage
#----------------------------#
if ($size_file) {
    my $length_of        = read_sizes($size_file);
    my $length_of_genome = sum0( values %{$length_of} );
    my @chrs             = sort keys %{$length_of};

    my @copies = sort { $a <=> $b } keys %{$count_of};
    my @low_copies = grep { $_ <= $low_cut and $_ > 1 } @copies;
    my %high_copies = map { $_ => 1 } grep { $_ > $low_cut } @copies;
    my $set_of = {};
    for my $copy (@copies) {
        $set_of->{$copy} = {};
        for my $chr (@chrs) {
            $set_of->{$copy}{$chr} = AlignDB::IntSpan->new;
            if ( exists $runlist_of->{$copy}{$chr} ) {
                $set_of->{$copy}{$chr}->add( $runlist_of->{$copy}{$chr} );
            }
        }
    }

    my $runlist2_of = { all => {}, };    # merged N and all
    my $sum_of      = { all => 0, };     # sum size of each copy
    my $piece_of    = { all => 0, };     # pieces of each copy
    for my $copy (@copies) {
        $runlist_of->{$copy} = {};
        if ( $high_copies{$copy} ) {     # there're changes that no high copies
            if ( !exists $runlist2_of->{N} ) {
                $runlist2_of->{N} = {};
            }
        }
        else {
            if ( $copy > 1 ) {
                $runlist2_of->{$copy} = {};
            }
        }

        for my $chr (@chrs) {
            $sum_of->{$copy} += $set_of->{$copy}{$chr}->count;

            # sum of all, there're overlaps between copies.
            $sum_of->{all}   += $sum_of->{$copy};
            $piece_of->{all} += $count_of->{$copy} * $copy;

            if ( !exists $runlist2_of->{all}{$chr} ) {
                $runlist2_of->{all}{$chr} = $set_of->{$copy}{$chr};
            }
            else {
                $runlist2_of->{all}{$chr}->add( $set_of->{$copy}{$chr} );
            }

            # sum of high copies
            if ( $high_copies{$copy} ) {
                $sum_of->{N}   += $sum_of->{$copy};
                $piece_of->{N} += $count_of->{$copy} * $copy;

                if ( !exists $runlist2_of->{N}{$chr} ) {
                    $runlist2_of->{N}{$chr} = $set_of->{$copy}{$chr};
                }
                else {
                    $runlist2_of->{N}{$chr}->add( $set_of->{$copy}{$chr} );
                }
            }
            else {
                if ( $copy > 1 ) {
                    $piece_of->{$copy} = $count_of->{$copy} * $copy;
                    $runlist2_of->{$copy}{$chr}
                        = $set_of->{$copy}{$chr};
                }
            }
        }
    }

    # for feather
    for my $i ( keys %{$runlist2_of} ) {
        for my $j ( keys %{ $runlist2_of->{$i} } ) {
            $runlist2_of->{$i}{$j} = $runlist2_of->{$i}{$j}->runlist;
        }
    }
    DumpFile( "$output.runlist.yml", $runlist2_of );

    my $coverage_of = {};
    for my $key ( keys %{$sum_of} ) {
        $coverage_of->{$key} = $sum_of->{$key} / $length_of_genome;
    }

    open my $fh, '>', "$output.csv";
    print {$fh} "copy,name,length,size,coverage,percent,count,piece,avg_size\n";
    for ( @low_copies, "all", ( exists $sum_of->{N} ? "N" : () ) ) {
        my $count;
        if ( $_ eq "all" ) {
            $count += $count_of->{$_} for @copies;
        }
        elsif ( $_ eq "N" ) {
            $count += $count_of->{$_} for keys %high_copies;
        }
        else {
            $count = $count_of->{$_};
        }
        print {$fh} join( ",",
            $_,                 "sum",
            $length_of_genome,  $sum_of->{$_},
            $coverage_of->{$_}, $sum_of->{$_} / $sum_of->{all},
            $count,             $piece_of->{$_},
            $sum_of->{$_} / $piece_of->{$_} ),
            "\n";
    }
    close $fh;
}

#----------------------------#
# write circos link file
#----------------------------#
{
    # linkN is actually hightlight file
    my $link_fh_of = {};
    for ( 2 .. $low_cut, 'N' ) {
        open my $fh, ">", "$output.link$_.txt";
        $link_fh_of->{$_} = $fh;
    }

    my @colors = reverse map {"paired-12-qual-$_"} ( 1 .. 12 );
    my $color_idx = 0;
    for my $c (@cc) {
        my $copy = scalar @{$c};
        next if $copy < 2;

        if ( $copy > $low_cut ) {
            for ( @{$c} ) {
                my ( $chr, $set, $strand ) = string_to_set($_);
                print { $link_fh_of->{N} } join( " ",
                    $chr, $set->min, $set->max,
                    "fill_color=" . $colors[$color_idx] ),
                    "\n";
            }

            # rotate color
            $color_idx++;
            $color_idx = 0 if $color_idx > 11;
            next;
        }

        for my $idx1 ( 0 .. $copy - 1 ) {
            for my $idx2 ( $idx1 + 1 .. $copy - 1 ) {
                my @fields;
                for ( $idx1, $idx2 ) {
                    my ( $chr, $set, $strand ) = string_to_set( $c->[$_] );
                    push @fields,
                        (
                        $chr,
                        $strand eq "+"
                        ? ( $set->min, $set->max )
                        : ( $set->max, $set->min )
                        );
                }
                print { $link_fh_of->{$copy} } join( " ", @fields ), "\n";
            }
        }
    }

    close $link_fh_of->{$_} for ( 2 .. $low_cut, 'N' );
}

#----------------------------#
# write piece sequences
#----------------------------#
if ( $size_file and $write_seq ) {
    my @chrs = sort keys %{ read_sizes($size_file) };

    my %file_of
        = map { $_ => file($size_file)->dir->file( $_ . ".fa" )->stringify }
        @chrs;
    my %seq_obj_of = map {
        $_ => Bio::SeqIO->new( -file => $file_of{$_}, -format => 'Fasta' )
            ->next_seq
    } @chrs;

    my $seq_fh_of = {};
    my @copies = grep { $_ >= 2 } sort { $a <=> $b } keys %{$count_of};
    for (@copies) {
        open my $fh, ">", "$output.copy$_.fas";
        $seq_fh_of->{$_} = $fh;
    }

    for my $c (@cc) {
        my $copy = scalar @{$c};
        next if $copy < 2;

        my @heads = @{$c};
        for my $head ( @{$c} ) {
            my ( $chr, $set, $strand ) = string_to_set($head);
            my $seq = $seq_obj_of{$chr}->subseq( $set->min, $set->max );
            $seq = uc $seq;
            if ( $strand ne "+" ) {
                $seq = revcom($seq);
            }

            print { $seq_fh_of->{$copy} } ">$head\n";
            print { $seq_fh_of->{$copy} } "$seq\n";
        }
        print { $seq_fh_of->{$copy} } "\n";
    }

    close $seq_fh_of->{$_} for (@copies);

}

$stopwatch->end_message;
exit;

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

sub string_to_set {
    my $node = shift;

    my ( $chr, $runlist ) = split /:/, $node;
    my $strand = "+";
    if ( $chr =~ /\((.+)\)/ ) {
        $strand = $1;
        $chr =~ s/\(.+\)//;
    }
    my $set = AlignDB::IntSpan->new($runlist);

    return ( $chr, $set, $strand );
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

    perl gather_info_axt.pl -f example.match.tsv

=cut


