#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use JSON::XS;
use File::Slurp;

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
    $output =~ s/\.lav$/\.json/ || die "Not a .lav file.\n";
}

#----------------------------------------------------------#
# now run!
#----------------------------------------------------------#
my $coder     = JSON::XS->new->ascii->pretty->allow_nonref;
my $full_json = {};

my $lav_content = read_file($lavfile);
my @lavs = split /\#\:lav/, $lav_content;
shift @lavs;    # .lav file start with #:lav

#----------------------------#
# d-stanza
#----------------------------#
{
    my $d_stanza = shift @lavs;
    $d_stanza =~ /d {\s+(.+?)\s+}/s;
    $d_stanza = $1;
    my @lines      = split /\n/, $d_stanza;
    my $line_first = shift @lines;
    my $line_last  = pop @lines;

    my %opt_of;
    my $opt_lines = $line_first . $line_last;
    while ( $opt_lines =~ /\b(\w)\s?=\s?(\d+)\b/g ) {
        $opt_of{$1} = $2;
    }
    $full_json->{header}{opt} = \%opt_of;

    my $score = join "\n", @lines;
    $full_json->{header}{score} = $score;
}

$full_json->{runs} = [];
foreach my $lav (@lavs) {
    my $run_json   = {};
    my $run_header = {};

    #----------------------------#
    # s-stanza
    #----------------------------#
    {
        $lav =~ /s {\s+(.+?)\s+}/s;
        my $s_stanza = $1;
        my @s_lines = $s_stanza =~ /(.+ \d+ \d+ \d+ \d+)/g;
        unless ( ( scalar @s_lines ) == 2 ) { die "s-stanza error.\n"; }

        $s_lines[0] =~ /\s*\"?(.+?)\-?\"? (\d+) (\d+) (\d+) (\d+)/;
        $run_header->{t_file}   = $1;
        $run_header->{t_from}   = $2;
        $run_header->{t_to}     = $3;
        $run_header->{t_strand} = $4 ? '-' : '+';
        $run_header->{t_number} = $5;

        $s_lines[1] =~ /\s*\"?(.+?)\-?\"? (\d+) (\d+) (\d+) (\d+)/;
        my ( $q_file, $q_strand, $q_contig ) = ( $1, $2, $3 );
        $run_header->{q_file}   = $1;
        $run_header->{q_from}   = $2;
        $run_header->{q_to}     = $3;
        $run_header->{q_strand} = $4 ? '-' : '+';
        $run_header->{q_number} = $5;
    }

    #----------------------------#
    # h-stanza
    #----------------------------#
    {
        $lav =~ /h {\s+(.+?)\s+}/s;
        my $h_stanza = $1;
        my @h_lines = $h_stanza =~ /(.+)/g;
        unless ( ( scalar @h_lines ) == 2 ) { die "h-stanza error.\n"; }

        $h_lines[0] =~ /\s*\"\>(.+)\"/;
        $run_header->{t_name} = $1;
        $run_header->{t_name} =~ s/ \(reverse complement\)//;

        $h_lines[1] =~ /\s*\"\>(.+)\"/;
        $run_header->{q_name} = $1;
        $run_header->{q_name} =~ s/ \(reverse complement\)//;
    }

    $run_json->{header} = $run_header;

    #----------------------------#
    # a-stanzas
    #----------------------------#
    $run_json->{aligns} = [];

    my @a_stanzas = $lav =~ /a {\s+(.+?)\s+}/sg;
    foreach my $a_stanza (@a_stanzas) {
        my $align_json = {};

        # s-line, scores
        unless ( $a_stanza =~ /\s*s (\d+)/ ) {
            die "There isn't a s-line.\n";
        }
        $align_json->{score} = $1;

        # b-line, begins
        unless ( $a_stanza =~ /\s*b (\d+) (\d+)/ ) {
            die "There isn't a s-line.\n";
        }
        $align_json->{t_from} = $1;
        $align_json->{q_from} = $2;

        # e-line, ends
        unless ( $a_stanza =~ /\s*e (\d+) (\d+)/ ) {
            die "There isn't a e-line.\n";
        }
        $align_json->{t_to} = $1;
        $align_json->{q_to} = $2;

        $align_json->{blocks} = [];

        my @align_pieces = $a_stanza =~ /\s*l (\d+ \d+ \d+ \d+) \d+/g;
        foreach my $align_piece (@align_pieces) {
            unless ( $align_piece =~ /(\d+) (\d+) (\d+) (\d+)/ ) {
                die "l-line error\n";
            }
            my $block_json = {};
            $block_json->{t_begin} = $1;
            $block_json->{q_begin} = $2;
            $block_json->{t_end}   = $3;
            $block_json->{q_end}   = $4;

            push @{ $align_json->{blocks} }, $block_json;
        }

        push @{ $run_json->{aligns} }, $align_json;
    }

    push @{ $full_json->{runs} }, $run_json;
}

my $text = $coder->encode($full_json);

open my $fh, '>', $output;
print {$fh} $text;
close $fh;

exit;

__END__

=head1 NAME

    lav2json.pl - convert .lav files to .json files

=head1 SYNOPSIS

    lav2json.pl -l <lavfile> -o <output>

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

