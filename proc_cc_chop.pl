#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;

use Bio::Seq;
use Bio::SeqIO;

use Set::Scalar;
use List::MoreUtils qw(minmax firstidx);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $cc_file;
my $size_file;

my $output;

my $msa = 'mafft';    # Default alignment program

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$cc_file,
    's|size=s'   => \$size_file,
    'o|output=s' => \$output,
    'msa=s'      => \$msa,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
if ( !defined $size_file ) {
    die "--size chr.sizes is needed\n";
}
elsif ( !-e $size_file ) {
    die "--size chr.sizes doesn't exist\n";
}

my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Analysis [$cc_file]");

if ( !$output ) {
    $output = path($cc_file)->basename;

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

my @cc       = @{ $yml->{cc} };
my $count_of = $yml->{count};
my @copies   = grep { $_ >= 2 } sort { $a <=> $b } keys %{$count_of};
    my @chrs = sort keys %{ read_sizes($size_file) };

#----------------------------#
# write piece sequences
#----------------------------#
{
    print "Write aligned cc sequences\n";

    print " " x 4, "Load genomic sequences\n";
    my %file_of
        = map { $_ => path($size_file)->parent->child( $_ . ".fa" )->stringify }
        @chrs;
    my %genome_of = map {
        $_ => Bio::SeqIO->new( -file => $file_of{$_}, -format => 'Fasta' )
            ->next_seq->seq
    } @chrs;

    # open file handlers
    my $seq_fh_of = {};
    for (@copies) {
        open my $fh, ">", "$output.copy$_.fas";
        $seq_fh_of->{$_} = $fh;
    }
    open $seq_fh_of->{pairwise}, ">", "$output.pairwise.fas";

    my @new_cc;

    print " " x 4, "Get cc sequences\n";
    for my $c (@cc) {
        my $copy = scalar @{$c};
        next if $copy < 2;
        print " " x 8, "This cc has $copy copies\n";

        # Get subseq from genomic files
        print " " x 12, "Get subseq from genomic files\n";
        my @heads = @{$c};
        my @seqs;
        for my $head ( @{$c} ) {
            my $seq = seq_from_string( \%genome_of, $head );
            push @seqs, $seq;
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

        my ( $head_chopped, $tail_chopped )
            = trim_head_tail( $seq_of, \@heads, 10, 10 );

        my ( $new_seq_of, $new_heads )
            = change_name_chopped( $seq_of, \@heads, $head_chopped,
            $tail_chopped );

        # names changed
        push @new_cc, $new_heads;

        # Write sequences
        print " " x 12, "Write sequences\n";
        for my $n ( @{$new_heads} ) {
            printf { $seq_fh_of->{$copy} } ">%s\n", $n;
            printf { $seq_fh_of->{$copy} } "%s\n",  $new_seq_of->{$n};
        }
        print { $seq_fh_of->{$copy} } "\n";

        # Write pairwise alignments
        print " " x 12, "Write pairwise alignments\n";
        my $pair_ary = best_pairwise( $new_seq_of, $new_heads );
        for my $p ( @{$pair_ary} ) {
            for my $n ( @{$p} ) {
                printf { $seq_fh_of->{pairwise} } ">%s\n", $n;
                printf { $seq_fh_of->{pairwise} } "%s\n",  $new_seq_of->{$n};
            }
            print { $seq_fh_of->{pairwise} } "\n";
        }
    }

    close $seq_fh_of->{$_} for (@copies);
    close $seq_fh_of->{pairwise};

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

DumpFile(
    "$output.chr.runlist.yml", $set_chr_of
);

$stopwatch->end_message;
exit;

sub seq_from_string {
    my $genome_of = shift;
    my $head      = shift;

    my ( $chr, $set, $strand ) = string_to_set($head);
    my $seq = substr $genome_of->{$chr}, $set->min - 1, $set->size;
    $seq = uc $seq;
    if ( $strand ne "+" ) {
        $seq = revcom($seq);
    }

    return $seq;
}

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
        my ( $min, $max ) = minmax @row;
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


