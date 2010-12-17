#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use File::Slurp;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $lavfile;
my $output;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'l|lavfile=s' => \$lavfile,
    'o|output=s'  => \$output,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$lavfile =~ s/\\/\//g;
unless ($output) {
    $output = $lavfile;
    $output =~ s/\.lav$/\.axt/ || die "Not a .lav file.\n";
}

#----------------------------------------------------------#
# now run!
#----------------------------------------------------------#
my $lav_content = read_file($lavfile);
my @lavs = split /\#\:lav/, $lav_content;
shift @lavs;    # .lav file start with #:lav
my $d_stanza = shift @lavs;    # Not needed by this program

open my $outfh, '>', $output;
my $align_id = 0;
my %cache;                     # cache fasta files
foreach my $lav (@lavs) {

    #----------------------------#
    # s-stanza
    #----------------------------#
    $lav =~ /s {\s+(.+?)\s+}/s;
    my $s_stanza = $1;
    my @s_lines = $s_stanza =~ /(.+ \s+ \d+ \s+ \d+ \s+ \d+ \s+ \d+)/gx;
    unless ( ( scalar @s_lines ) == 2 ) { die "s-stanza error.\n"; }

    $s_lines[0] =~ /\s*\"?(.+?)\-?\"? \s+ \d+ \s+ \d+ \s+ (\d+) \s+ (\d+)/x;
    my ( $t_file, $t_strand, $t_contig ) = ( $1, $2, $3 );

    $s_lines[1] =~ /\s*\"?(.+?)\-?\"? \s+ \d+ \s+ \d+ \s+ (\d+) \s+ (\d+)/x;
    my ( $q_file, $q_strand, $q_contig ) = ( $1, $2, $3 );

    #----------------------------#
    # h-stanza
    #----------------------------#
    $lav =~ /h {\s+(.+?)\s+}/s;
    my $h_stanza = $1;
    my @h_lines = $h_stanza =~ /(.+)/g;
    unless ( ( scalar @h_lines ) == 2 ) { die "h-stanza error.\n"; }

    #----------------------------#
    # generate two sequences
    #----------------------------#
    if ( !exists $cache{$t_file} ) {
        my ( $seq_of, $seq_names ) = read_fasta($t_file);
        $cache{$t_file} = {
            seq_of    => $seq_of,
            seq_names => $seq_names,
        };
    }
    my $t_name = $cache{$t_file}->{seq_names}->[ $t_contig - 1 ];
    my $t_seq  = $cache{$t_file}->{seq_of}->{$t_name};

    if ( !exists $cache{$q_file} ) {
        my ( $seq_of, $seq_names ) = read_fasta($q_file);
        $cache{$q_file} = {
            seq_of    => $seq_of,
            seq_names => $seq_names,
        };
    }
    my $q_name = $cache{$q_file}->{seq_names}->[ $q_contig - 1 ];
    my $q_seq  = $cache{$q_file}->{seq_of}->{$q_name};
    if ($q_strand) {
        $q_seq = revcom($q_seq);
    }

    #----------------------------#
    # generate axt alignments
    #----------------------------#
    my @a_stanzas = $lav =~ /a {\s+(.+?)\s+}/sg;
    foreach my $a_stanza (@a_stanzas) {
        my $alignment_target  = '';
        my $alignment_query   = '';
        my @align_pieces      = $a_stanza =~ /\s*l (\d+ \d+ \d+ \d+) \d+/g;
        my $former_end_target = '';
        my $former_end_query  = '';
        foreach my $align_piece (@align_pieces) {
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
            $alignment_target .= (
                substr $t_seq,
                ( $t_begin - 1 - ( length $q_del ) ),
                $length_target
            );
            $alignment_query .= (
                substr $q_seq,
                ( $q_begin - 1 - ( length $t_del ) ),
                $length_query
            );
            if ( ( length $alignment_query ) ne ( length $alignment_target ) )
            {
                die "Target length doesn't match query's in the alignment.\n";
            }
            $former_end_target = $t_end;
            $former_end_query  = $q_end;
        }

        # b-line, begins
        unless ( $a_stanza =~ /\s*b (\d+) (\d+)/ ) {
            die "There isn't a s-line.\n";
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

=head1 NAME

    lav2axt.pl - convert .lav files to .axt files

=head1 SYNOPSIS

    lav2axt.pl -l <lavfile> -o <output>

    lav2axt.pl [options]
     Options:
       -h, --help               brief help message
       -m, --man                full documentation
       -l, --lavfile
       -o, --output

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

