#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;

use Bio::SeqIO;

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

#----------------------------#
# write piece sequences
#----------------------------#
{
    print "Write aligned cc sequences\n";
    my @chrs = sort keys %{ read_sizes($size_file) };

    print " " x 4, "Load genomic sequences\n";
    my %file_of
        = map { $_ => path($size_file)->parent->child( $_ . ".fa" )->stringify }
        @chrs;
    my %genome_of = map {
        $_ => Bio::SeqIO->new( -file => $file_of{$_}, -format => 'Fasta' )
            ->next_seq->seq
    } @chrs;

    my $seq_fh_of = {};
    for (@copies) {
        open my $fh, ">", "$output.copy$_.fas";
        $seq_fh_of->{$_} = $fh;
    }

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

        my ( $neq_seq_of, $new_heads )
            = change_name_chopped( $seq_of, \@heads, $head_chopped,
            $tail_chopped );

        # names changed
        push @new_cc, $new_heads;

        # Write sequences
        print " " x 12, "Write sequences\n";
        for my $n ( @{$new_heads} ) {
            printf { $seq_fh_of->{$copy} } ">%s\n", $n;
            printf { $seq_fh_of->{$copy} } "%s\n",  $neq_seq_of->{$n};
        }
        print { $seq_fh_of->{$copy} } "\n";
    }

    close $seq_fh_of->{$_} for (@copies);

    # we don't sort @new_cc
    @cc = @new_cc;
    print "\n";
}

#----------------------------#
# runlist
#----------------------------#
my $set_of = {};
for my $c (@cc) {
    my $copy = scalar @{$c};
    if ( !exists $set_of->{$copy} ) {
        $set_of->{$copy} = {};
    }
    for ( @{$c} ) {
        my ( $chr, $set, $strand ) = string_to_set($_);
        if ( !exists $set_of->{$copy}{$chr} ) {
            $set_of->{$copy}{$chr} = AlignDB::IntSpan->new;
        }
        $set_of->{$copy}{$chr}->add($set);
    }
}

for my $key_i ( keys %{$set_of} ) {
    for my $key_j ( keys %{ $set_of->{$key_i} } ) {
        $set_of->{$key_i}{$key_j} = $set_of->{$key_i}{$key_j}->runlist;
    }
}

print "Write chopped cc\n";
DumpFile(
    "$output.yml",
    {   count   => $count_of,
        cc      => \@cc,
        runlist => $set_of,
    }
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


