#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use App::Fasops::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

lav2axt.pl - convert .lav files to .axt files

=head1 SYNOPSIS

    perl lav2axt.pl -i <lavfile(s)> -o <output> [options]
      Options:
        --help          -?          brief help message
        --input         -i  STR     input lav file
        --output        -o  STR     output axt file

=head1 lav format

    http://www.bx.psu.edu/miller_lab/dist/lav_format.html

Here <start> and <stop> are origin 1 (i.e. the first base in the original given sequence is called
"1") and inclusive (both endpoints are included in the interval).

=head1 axt format

    https://genome.ucsc.edu/goldenPath/help/axt.html

If the strand value is "-", the values of the aligning organism's start and end fields are relative
to the reverse-complemented coordinates of its chromosome.

=cut

GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'input|i=s'  => \my $lavfile,
    'output|o=s' => \my $output,
) or Getopt::Long::HelpMessage(1);

$lavfile =~ s/\\/\//g;
unless ($output) {
    $output = $lavfile;
    $output =~ s/\.lav$/\.axt/ or die "Not a .lav file.\n";
}

#----------------------------------------------------------#
# now run!
#----------------------------------------------------------#
my $lav_content = path($lavfile)->slurp;
my @lavs = split /\#\:lav/, $lav_content;
shift @lavs;    # .lav file start with #:lav
my $d_stanza = shift @lavs;    # Not needed by this program

open my $outfh, '>', $output;
my $align_id = 0;
my %cache;                     # cache fasta files
for my $lav (@lavs) {

    #----------------------------#
    # s-stanza
    #----------------------------#
    # "<filename>[-]" <start> <stop> [<rev_comp_flag> <sequence_number>]
    $lav =~ /s {\s+(.+?)\s+}/s;
    my $s_stanza = $1;
    my @s_lines  = $s_stanza =~ /(.+ \s+ \d+ \s+ \d+ \s+ \d+ \s+ \d+)/gx;
    if ( scalar @s_lines != 2 ) {
        die "s-stanza error.\n";
    }

    $s_lines[0] =~ /\s*\"?(.+?)\-?\"? \s+ (\d+) \s+ (\d+) \s+ (\d+) \s+ (\d+)/x;
    my ( $t_file, $t_seq_start, $t_seq_stop, $t_strand, $t_contig )
        = ( $1, $2, $3, $4, $5 );
    if ( $t_seq_start != 1 ) {
        die "Target sequence doesn't start at 1\n";
    }

    $s_lines[1] =~ /\s*\"?(.+?)\-?\"? \s+ (\d+) \s+ (\d+) \s+ (\d+) \s+ (\d+)/x;
    my ( $q_file, $q_seq_start, $q_seq_stop, $q_strand, $q_contig )
        = ( $1, $2, $3, $4, $5 );
    if ( $q_seq_start != 1 ) {
        die "Query sequence doesn't start at 1\n";
    }

    #----------------------------#
    # h-stanza
    #----------------------------#
    $lav =~ /h {\s+(.+?)\s+}/s;
    my $h_stanza = $1;
    my @h_lines  = $h_stanza =~ /(.+)/g;
    if ( scalar @h_lines != 2 ) {
        die "h-stanza error.\n";
    }

    #----------------------------#
    # generate two sequences
    #----------------------------#
    if ( !exists $cache{$t_file} ) {
        my $seq_of = App::Fasops::Common::read_fasta($t_file);
        $cache{$t_file} = {
            seq_of    => $seq_of,
            seq_names => keys %{$seq_of},
        };
    }
    my $t_name = $cache{$t_file}->{seq_names}->[ $t_contig - 1 ];
    my $t_seq  = $cache{$t_file}->{seq_of}->{$t_name};

    if ( !exists $cache{$q_file} ) {
        my $seq_of = App::Fasops::Common::read_fasta($q_file);
        $cache{$q_file} = {
            seq_of    => $seq_of,
            seq_names => keys %{$seq_of},
        };
    }
    my $q_name = $cache{$q_file}->{seq_names}->[ $q_contig - 1 ];
    my $q_seq  = $cache{$q_file}->{seq_of}->{$q_name};
    if ($q_strand) {
        $q_seq = App::Fasops::Common::revcom($q_seq);
    }

    #----------------------------#
    # generate axt alignments
    #----------------------------#
    my @a_stanzas = $lav =~ /a {\s+(.+?)\s+}/sg;
    for my $a_stanza (@a_stanzas) {
        my $alignment_target  = '';
        my $alignment_query   = '';
        my @align_pieces      = $a_stanza =~ /\s*l (\d+ \d+ \d+ \d+) \d+/g;
        my $former_end_target = '';
        my $former_end_query  = '';
        for my $align_piece (@align_pieces) {
            unless ( $align_piece =~ /(\d+) (\d+) (\d+) (\d+)/g ) {
                die "l-line error\n";
            }
            my ( $t_begin, $q_begin, $t_end, $q_end ) = ( $1, $2, $3, $4 );
            my $t_del = '';
            my $q_del = '';
            if ( $alignment_target
                && ( $t_begin - $former_end_target > 1 ) )
            {
                $q_del = '-' x ( $t_begin - $former_end_target - 1 );
                $alignment_query .= $q_del;
            }
            if ( $alignment_query
                && ( $q_begin - $former_end_query > 1 ) )
            {
                $t_del = '-' x ( $q_begin - $former_end_query - 1 );
                $alignment_target .= $t_del;
            }
            my $length_target = $t_end - $t_begin + 1 + ( length $q_del );
            my $length_query  = $q_end - $q_begin + 1 + ( length $t_del );
            $alignment_target
                .= substr( $t_seq, ( $t_begin - 1 - ( length $q_del ) ),
                $length_target );
            $alignment_query
                .= substr( $q_seq, ( $q_begin - 1 - ( length $t_del ) ),
                $length_query );
            if ( ( length $alignment_query ) ne ( length $alignment_target ) ) {
                die "Target length doesn't match query's in the alignment.\n";
            }
            $former_end_target = $t_end;
            $former_end_query  = $q_end;
        }

        # b-line, begins
        unless ( $a_stanza =~ /\s*b (\d+) (\d+)/ ) {
            die "There isn't a b-line.\n";
        }
        my $t_from = $1;
        my $q_from = $2;

        # e-line, ends
        unless ( $a_stanza =~ /\s*e (\d+) (\d+)/ ) {
            die "There isn't a e-line.\n";
        }
        my $t_to = $1;
        my $q_to = $2;

        # s-line, scores
        unless ( $a_stanza =~ /\s*s (\d+)/ ) {
            die "There isn't a s-line.\n";
        }
        my $score = $1;

        # only keep the first part in fasta header
        ($t_name) = split /\W+/, $t_name;
        ($q_name) = split /\W+/, $q_name;

        # prepare axt header
        my $axt_head = $align_id;
        $axt_head .= " $t_name $t_from $t_to";
        $axt_head .= " $q_name $q_from $q_to ";
        $axt_head .= $q_strand ? '-' : '+';
        $axt_head .= " $score\n";

        # write axt file
        print {$outfh} $axt_head;
        print {$outfh} "$alignment_target\n";
        print {$outfh} "$alignment_query\n\n";
        $align_id++;
    }
}

close $outfh;
exit;

__END__
