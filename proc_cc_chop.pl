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

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(multi_align multi_align_matrix trim_head_tail);

use lib "$FindBin::RealBin/lib";
use MyUtil qw(string_to_set revcom get_seq_faidx read_sizes decode_header change_name_chopped);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

merge_node.pl - merge overlapped nodes of paralog graph
    
=head1 SYNOPSIS

    perl merge_node.pl -f <file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     file
        --size          -s  STR     chr.sizes
        --genome        -g  STR     reference genome file
        --output        -o  STR     output   
        --msa               STR     Aligning program, default is [mafft]

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \my $cc_file,
    'size|s=s'   => \my $size_file,
    'genome|g=s' => \my $genome,
    'output|o=s' => \my $output,
    'msa=s' => \( my $msa = 'mafft' ),
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
    $output = "$output.cc";
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
my $yml = LoadFile($cc_file);

my @cc       = @{ $yml->{cc} };
my $count_of = $yml->{count};
my @copies   = grep { $_ >= 2 } sort { $a <=> $b } keys %{$count_of};
my @chrs     = keys %{ read_sizes($size_file) };

#----------------------------#
# write piece sequences
#----------------------------#
{
    print "Write aligned cc sequences\n";

    for (@copies) {
        path("$output.copy$_.fas")->remove;
    }
    path("$output.pairwise.fas")->remove;

    my @new_cc;
    
    my @bad_cc;

    print " " x 4 . "Get cc sequences\n";
    for my $c (@cc) {
        my $copy = scalar @{$c};
        if ( $copy < 2 ) {
            printf " " x 8 . "Single cc [%s]\n", $c->[0];
            push @bad_cc, $c;
            next;
        }
        print " " x 8 . "This cc has $copy copies\n";

        # Get subseq from genomic files
        print " " x 12 . "Get sequences from genomic files\n";
        my @heads = @{$c};
        my @seqs;
        for my $head ( @{$c} ) {
            my $info = decode_header($head);
            my $location = sprintf "%s:%d-%d", $info->{chr_name}, $info->{chr_start},
                $info->{chr_end};
            my $seq = get_seq_faidx( $genome, $location );
            printf " " x 16 . "Location: %s; Length %d\n", $location, length $seq;
            if ( $info->{chr_strand} ne "+" ) {
                $seq = revcom($seq);
            }
            push @seqs, $seq;
        }
        my ( $l_min, $l_max ) = minmax( map { length $_ } @seqs );
        my $diff_ratio = sprintf "%.3f", $l_min / $l_max;
        print " " x 12 . "Ratio[$diff_ratio]\tMin [$l_min]\tMax[$l_max]\n";

        if ( $diff_ratio < 0.8 ) {
            print " " x 12 . "*** Don't align this cc. Check it.\n";
            push @bad_cc, $c;
            next;
        }

        # aligning
        print " " x 12, "Realign with $msa\n";
        my $realigned_seqs = multi_align( \@seqs, $msa );

        # Chopping hanging parts
        print " " x 12, "Chopping hanging parts\n";
        my $seq_of = {};
        for my $i ( 0 .. $#heads ) {
            $seq_of->{ $heads[$i] } = uc $realigned_seqs->[$i];
        }

        my ( $head_chopped, $tail_chopped ) = trim_head_tail( $seq_of, \@heads, 10, 10 );

        my ( $new_seq_of, $new_heads )
            = change_name_chopped( $seq_of, \@heads, $head_chopped, $tail_chopped );
        printf " " x 16 . "Chopping %d %d\n", max( values %$head_chopped ),
            max( values %$tail_chopped );

        # names changed
        push @new_cc, $new_heads;

        # Write sequences
        print " " x 12, "Write sequences\n";
        for my $n ( @{$new_heads} ) {
            path("$output.copy$copy.fas")->append(">$n\n");
            path("$output.copy$copy.fas")->append( $new_seq_of->{$n} . "\n" );
        }
        path("$output.copy$copy.fas")->append("\n");

        # Write pairwise alignments
        print " " x 12, "Write pairwise alignments\n";
        my $pair_ary = best_pairwise( $new_seq_of, $new_heads );
        for my $p ( @{$pair_ary} ) {
            my $realigned_pair
                = multi_align( [ $new_seq_of->{ $p->[0] }, $new_seq_of->{ $p->[1] } ], $msa );
            for my $i ( 0 .. scalar @{$p} - 1 ) {
                path("$output.pairwise.fas")->append( ">" . $p->[$i] . "|copy=" . $copy . "\n" );
                path("$output.pairwise.fas")->append( uc( $realigned_pair->[$i] ) . "\n" );
            }
            path("$output.pairwise.fas")->append("\n");
        }
    }

    DumpFile("BAD_cc.yml", \@bad_cc);
    
    # we don't sort @new_cc
    @cc = @new_cc;
    print "\n";
    
}

#----------------------------#
# runlist
#----------------------------#
# per copy and chr
my $set_copy_of = {};
for my $c (@cc) {
    my $copy = scalar @{$c};
    if ( !exists $set_copy_of->{$copy} ) {
        $set_copy_of->{$copy} = {};
    }
    for ( @{$c} ) {
        my ( $chr, $set, $strand ) = string_to_set($_);
        if ( !exists $set_copy_of->{$copy}{$chr} ) {
            $set_copy_of->{$copy}{$chr} = AlignDB::IntSpan->new;
        }
        $set_copy_of->{$copy}{$chr}->add($set);
    }
}

for my $key_i ( keys %{$set_copy_of} ) {
    for my $key_j ( keys %{ $set_copy_of->{$key_i} } ) {
        $set_copy_of->{$key_i}{$key_j}
            = $set_copy_of->{$key_i}{$key_j}->runlist;
    }
}

print "Write chopped cc\n";
DumpFile(
    "$output.yml",
    {   count   => $count_of,
        cc      => \@cc,
        runlist => $set_copy_of,
    }
);

# per chr
my $set_chr_of = {};

for my $c (@cc) {
    for ( @{$c} ) {
        my ( $chr, $set, $strand ) = string_to_set($_);
        if ( !exists $set_chr_of->{$chr} ) {
            $set_chr_of->{$chr} = AlignDB::IntSpan->new;
        }
        $set_chr_of->{$chr}->add($set);
    }
}
for my $key_i ( keys %{$set_chr_of} ) {
    $set_chr_of->{$key_i} = $set_chr_of->{$key_i}->runlist;
}

DumpFile( "$output.chr.runlist.yml", $set_chr_of );

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

__END__
