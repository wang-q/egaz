#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use Set::Scalar;
use List::Util qw(max);
use List::MoreUtils qw(minmax firstidx);

use MCE;
use MCE::Flow Sereal => 1;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(multi_align multi_align_matrix);

use lib "$FindBin::RealBin/lib";
use MyUtil
    qw(string_to_set revcom get_seq_faidx read_sizes decode_header change_name_chopped sort_cc);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

proc_cc_chop.pl - chopping connected paralog sequences
    
=head1 SYNOPSIS

    perl proc_cc_chop.pl -f <cc file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     file
        --size          -s  STR     chr.sizes
        --genome        -g  STR     reference genome file
        --ratio         -r  FLOAT   report bad cc if lengths differences
                                    large than this ratio, default is [0.8]
        --output        -o  STR     output
        --msa               STR     Aligning program, default is [mafft]
        --parallel      -p  INT     default is [8]

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \( my $cc_file ),
    'size|s=s'   => \( my $size_file ),
    'genome|g=s' => \( my $genome ),
    'ratio|r=f' => \( my $ratio = 0.8 ),
    'output|o=s'   => \( my $output ),
    'msa=s'        => \( my $msa = 'mafft' ),
    'parallel|p=i' => \( my $parallel = 8 ),
) or HelpMessage(1);

if ( !defined $size_file ) {
    die "--size is needed\n";
}
elsif ( !path($size_file)->is_file ) {
    die "--size doesn't exist\n";
}

if ( !defined $genome ) {
    die "--genome is needed\n";
}
elsif ( !path($genome)->is_file ) {
    die "--genome doesn't exist\n";
}

if ( !$output ) {
    $output = path($cc_file)->basename;
    ($output) = grep {defined} split /\./, $output;
}

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$cc_file]");

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

#----------------------------#
# load cc
#----------------------------#
$stopwatch->block_message("Load cc [$cc_file]");
my @cc           = @{ LoadFile($cc_file) };
my $total_number = scalar @cc;

{    # Remove alignment files if existing
    my %seen;
    for my $c (@cc) {
        my $copy = scalar @{$c};
        $seen{$copy}++;
    }

    for ( keys %seen ) {
        path("$output.copy$_.fas")->remove;
    }
    path("$output.pairwise.fas")->remove;
}

#----------------------------#
# write piece sequences
#----------------------------#
$stopwatch->block_message("Write aligned cc sequences");

my $file_bad_cc = "bad_cc.txt";
path($file_bad_cc)->remove;

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $c    = $chunk_ref->[0];
    my $wid  = MCE->wid;
    my $copy = scalar @{$c};

    print "* Process cc [$chunk_id/$total_number] with [$copy] copies by worker #$wid\n";

    # Get subseq from genomic files
    print " " x 4 . "Get sequences from genomic files\n";
    my @heads = @{$c};
    my @seqs;
    for my $head ( @{$c} ) {
        my $info     = decode_header($head);
        my $location = sprintf "%s:%d-%d", $info->{chr_name}, $info->{chr_start}, $info->{chr_end};
        my $seq      = get_seq_faidx( $genome, $location );
        printf " " x 8 . "Location: %s; Length %d\n", $location, length $seq;
        if ( $info->{chr_strand} ne "+" ) {
            $seq = revcom($seq);
        }
        push @seqs, $seq;
    }
    my ( $l_min, $l_max ) = minmax( map { length $_ } @seqs );
    my $diff_ratio = sprintf "%.3f", $l_min / $l_max;
    print " " x 4 . "Ratio[$diff_ratio]\tMin [$l_min]\tMax[$l_max]\n";

    if ( $diff_ratio < $ratio ) {
        print " " x 4 . "*** Don't align this cc. Check it.\n";
        path($file_bad_cc)->append( Dump($c) );
        next;
    }

    # aligning
    print " " x 4, "Realign with $msa\n";
    my $realigned_seqs = multi_align( \@seqs, $msa );

    # Chopping hanging parts
    print " " x 4, "Chopping hanging parts\n";
    my $seq_of = {};
    for my $i ( 0 .. $#heads ) {
        $seq_of->{ $heads[$i] } = uc $realigned_seqs->[$i];
    }

    my ( $head_chopped, $tail_chopped ) = trim_head_tail( $seq_of, \@heads, 10, 10 );

    printf " " x 8 . "Chopping %d %d\n", max( values %$head_chopped ), max( values %$tail_chopped );
    my ( $new_seq_of, $new_heads )
        = change_name_chopped( $seq_of, \@heads, $head_chopped, $tail_chopped );

    # Write sequences
    print " " x 4, "Write sequences\n";
    my $str = "";
    for my $n ( @{$new_heads} ) {
        $str = ">$n\n";
        $str .= uc( $new_seq_of->{$n} ) . "\n";
    }
    $str .= "\n";
    path("$output.cc.copy$copy.fas")->append($str);

    # Write pairwise alignments
    print " " x 4, "Write pairwise alignments\n";
    my $pair_ary = best_pairwise( $new_seq_of, $new_heads );
    $str = "";
    for my $p ( @{$pair_ary} ) {
        my $realigned_pair
            = multi_align( [ $new_seq_of->{ $p->[0] }, $new_seq_of->{ $p->[1] } ], $msa );
        for my $i ( 0 .. scalar @{$p} - 1 ) {
            $str .= ">" . $p->[$i] . "|copy=" . $copy . "\n";
            $str .= uc( $realigned_pair->[$i] ) . "\n";
        }
        $str .= "\n";
    }
    path("$output.cc.pairwise.fas")->append($str);

    # names changed
    print " " x 4, "Gather new heads for cc\n";
    MCE->gather($new_heads);
};

MCE::Flow::init {
    chunk_size  => 1,
    max_workers => $parallel,
};
my @cc_new = mce_flow $worker, \@cc;
MCE::Flow::finish;

$stopwatch->block_message("Write file [$output.cc.yml]");
@cc_new = sort_cc(@cc_new);
DumpFile( "$output.cc.yml", \@cc_new );

if ( path($file_bad_cc)->is_file ) {
    $stopwatch->block_message("There are bad cc. Check [$file_bad_cc]");
}

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub best_pairwise {
    my $seq_of    = shift;
    my $seq_names = shift;

    my $seq_number = scalar @{$seq_names};

    my @seqs = map { $seq_of->{$_} } @{$seq_names};

    my $matrix = multi_align_matrix( \@seqs );
    my $values = $matrix->_values;

    my $pair_set = Set::Scalar->new;
    for ( my $i = 0; $i < $seq_number; $i++ ) {
        my @row = @{ $values->[$i] };
        splice @row, $i, 1;    # remove the score of this item
        my ( $min, $max ) = minmax(@row);
        my $min_idx = firstidx { $_ == $min } @{ $values->[$i] };
        my @pair = ( $i, $min_idx );
        @pair = sort { $a <=> $b } @pair;
        $pair_set->insert( $pair[0] . ':' . $pair[1] );
    }

    my @pair_ary;
    for my $m ( $pair_set->members ) {
        my @pair = split ':', $m;
        push @pair_ary, [ $seq_names->[ $pair[0] ], $seq_names->[ $pair[1] ] ];
    }

    return \@pair_ary;
}


#----------------------------#
# trim head and tail indels
#----------------------------#
#  If head length set to 1, the first indel will be trimmed
#  Length set to 5 and the second indel will also be trimmed
#   GAAA--C
#   --AAAGC
#   GAAAAGC
sub trim_head_tail {
    my $seq_of      = shift;
    my $seq_names   = shift;
    my $head_length = shift;    # indels in this region will also be trimmed
    my $tail_length = shift;    # indels in this region will also be trimmed

    # default value means only trim indels starting at the first base
    $head_length = defined $head_length ? $head_length : 1;
    $tail_length = defined $tail_length ? $tail_length : 1;

    my $seq_number   = scalar @{$seq_names};
    my $align_length = length $seq_of->{ $seq_names->[0] };

    my $align_set = AlignDB::IntSpan->new("1-$align_length");
    my $indel_set = AlignDB::IntSpan->new;

    for my $n ( @{$seq_names} ) {
        my $seq_indel_set = find_indel_set( $seq_of->{$n} );
        $indel_set->merge($seq_indel_set);
    }

    # record bp chopped
    my %head_chopped = map { $_ => 0 } @{$seq_names};
    my %tail_chopped = map { $_ => 0 } @{$seq_names};

    # There're no indels at all
    return ( \%head_chopped, \%tail_chopped ) if $indel_set->is_empty;

    # head indel(s) to be trimmed
    my $head_set = AlignDB::IntSpan->new;
    $head_set->add_range( 1, $head_length );
    my $head_indel_set = $indel_set->find_islands($head_set);

    # head indels
    if ( $head_indel_set->is_not_empty ) {
        for my $i ( 1 .. $head_indel_set->max ) {
            my @column;
            for my $n ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$n}, 0, 1, '' );
                if ( $base ne '-' ) {
                    $head_chopped{$n}++;
                }
            }
        }
    }

    # tail indel(s) to be trimmed
    my $tail_set = AlignDB::IntSpan->new;
    $tail_set->add_range( $align_length - $tail_length + 1, $align_length );
    my $tail_indel_set = $indel_set->find_islands($tail_set);

    # tail indels
    if ( $tail_indel_set->is_not_empty ) {
        for my $i ( $tail_indel_set->min .. $align_length ) {
            my @column;
            for my $n ( @{$seq_names} ) {
                my $base = substr( $seq_of->{$n}, -1, 1, '' );
                if ( $base ne '-' ) {
                    $tail_chopped{$n}++;
                }
            }
        }
    }

    return ( \%head_chopped, \%tail_chopped );
}

__END__
