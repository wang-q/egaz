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

use FindBin;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $file;
my $size_file;

my $output;
my $outdir;

my $low_cut = 4;    # low copies their own stats

my $verbose;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$file,
    's|size=s'   => \$size_file,
    'o|output=s' => \$output,
    'd|outdir=s' => \$outdir,
    'l|low=s'    => \$low_cut,
    'v|verbose'  => \$verbose,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$file]");

if ( !$output ) {
    $output = basename($file);

    ($output) = grep {defined} split /\./, $output;
    $output = "$output.cc";
}
$output = File::Spec->catfile( $outdir, $output ) if $outdir;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $g = LoadFile($file);

#----------------------------#
# create cc
#----------------------------#
my @cc = $g->connected_components;

# sort by chromosome order within cc
for my $c (@cc) {
    my @unsorted = @{$c};

    # start point on chromosomes
    @unsorted = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @unsorted;

    # chromosome name
    @unsorted = map { $_->[0] }
        sort { $a->[1] cmp $b->[1] }
        map { /([\w.]+)\(.\)\:/; [ $_, $1 ] } @unsorted;

    $c = [@unsorted];
}

# sort by first node's chromosome order between cc
@cc = map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    map { $_->[0] =~ /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @cc;

@cc = map { $_->[0] }
    sort { $a->[1] cmp $b->[1] }
    map { $_->[0] =~ /([\w.]+)\(.\)\:/; [ $_, $1 ] } @cc;

# sort by nodes number between cc
@cc = sort { scalar @{$b} <=> scalar @{$a} } @cc;

#----------------------------#
# count
#----------------------------#
my $count_of = {};
for my $c (@cc) {
    my $copy = scalar @{$c};
    $count_of->{$copy}++;
}
DumpFile( "$output.yml", { count => $count_of, cc => \@cc, } );

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
        }
    }
    for my $c (@cc) {
        my $copy = scalar @{$c};
        for ( @{$c} ) {
            my ( $chr, $set, $strand ) = string_to_set($_);
            $set_of->{$copy}{$chr}->add($set);
        }
    }

    my $runlist_of  = {};    # original runlists of each copy
    my $runlist2_of = {};    # merge N
    my $sum_of      = {};    # sum size of each copy
    my $piece_of    = {};    # pieces of each copy
    for my $copy (@copies) {
        $runlist_of->{$copy} = {};
        if ( $high_copies{$copy} ) {
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
            $runlist_of->{$copy}{$chr} = $set_of->{$copy}{$chr}->runlist;
            $sum_of->{$copy} += $set_of->{$copy}{$chr}->count;

            # sum of all, there're overlaps between copies.
            $sum_of->{all} += $sum_of->{$copy};

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

    # override "$output.yml"
    DumpFile( "$output.yml",
        { count => $count_of, cc => \@cc, runlist => $runlist_of, } );

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
    for ( @low_copies, ( exists $sum_of->{N} ? "N" : () ) ) {
        my $count;
        if ( $_ eq "N" ) {
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
# linkN is actually hightlight file
my $fh_of = {};
for ( 2 .. $low_cut, 'N' ) {
    open my $fh, ">", "$output.link$_.txt";
    $fh_of->{$_} = $fh;
}

my @colors = reverse map {"paired-12-qual-$_"} ( 1 .. 12 );
my $color_idx = 0;
for my $c (@cc) {
    my $copy = scalar @{$c};
    next if $copy < 2;

    if ( $copy > $low_cut ) {
        for ( @{$c} ) {
            my ( $chr, $set, $strand ) = string_to_set($_);
            print { $fh_of->{N} } join( " ",
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
            next unless $g->has_edge( $c->[$idx1], $c->[$idx2] );
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
            print { $fh_of->{$copy} } join( " ", @fields ), "\n";
        }
    }
}

close $fh_of->{$_} for ( 2 .. $low_cut, 'N' );

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


