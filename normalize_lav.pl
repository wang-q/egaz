#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use List::Util qw(max);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'input|i=s'  => \my $lavfile,
    'output|o=s' => \my $outfile,
    'len0|0=i'   => \(my $len0 = 0),
    'len1|1=i'   => \(my $len1 = 0),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# now run!
#----------------------------------------------------------#
open my $fh, '>', $outfile;

my $lav_content = path($lavfile)->slurp;
my @lavs = split /\#\:lav\n/, $lav_content;
shift @lavs;    # .lav file start with #:lav
my $d_stanza = shift @lavs;
$d_stanza = "d {\n  normalize-lav $len0 $len1\n}\n" . $d_stanza;
print {$fh} "#:lav\n";
print {$fh} $d_stanza;

for my $lav (@lavs) {
    print {$fh} "#:lav\n";

    my $t_from = 0;
    my $q_from = 0;
    my $t_to   = 0;
    my $q_to   = 0;
    my $isrc   = 0;

    #----------------------------#
    # s-stanza
    #----------------------------#
    $lav =~ /s {\s+(.+?)\s+}/s;
    my $s_stanza = $1;
    my @s_lines = $s_stanza =~ /(.+ \d+ \d+ \d+ \d+)/g;
    unless ( ( scalar @s_lines ) == 2 ) { die "s-stanza error.\n"; }

    print {$fh} "s {\n";

    $s_lines[0] =~ /^\s*("[^"]*")\s+(\d+)\s+(\d+)\s+(.*)$/ or die;
    $t_from = $2;
    $t_to   = $3;
    print {$fh} "  $1 1 " . max( $t_to, $len0 ) . " $4\n";

    $s_lines[1] =~ /^\s*("[^"]*")\s+(\d+)\s+(\d+)\s+(.*)$/ or die;
    $q_from = $2;
    $q_to   = $3;
    print {$fh} "  $1 1 " . max( $q_to, $len1 ) . " $4\n";

    $isrc = scalar( $1 =~ /-"$/ );

    print {$fh} "}\n";

    #----------------------------#
    # h-stanza
    #----------------------------#
    if ( $lav =~ /h {\n(.+?)}/s ) {
        my $h_stanza = $1;
        print {$fh} "h {\n$h_stanza}\n";
    }

    #----------------------------#
    # a-stanza
    #----------------------------#
    # abs: 1....from....to....len1
    # rel: .....1.......n.........
    #
    # n == (to-from+1)
    # rev(x) == n-x+1
    # abs(x) == from+x-1
    #
    # abs(rev(x)) == from+(to-from+1-x+1)-1 == to-x+1
    # Rev(abs(rev(x))) == len1-(to-x+1)+1 == len1-to+x

    my $S0 = sub { $_[0] + $t_from - 1; };
    my $M1 = sub {
        my $x = $_[0];
        if ( $len1 == 0 ) {
            return $x;
        }
        else {
            return $isrc ? $len1 - $q_to + $x : $q_from + $x - 1;
        }
    };

    my @a_stanzas = $lav =~ /a {\s+(.+?)\s+}/sg;
    for my $a_stanza (@a_stanzas) {
        print {$fh} "a {\n";
        for my $line ( split /\n/, $a_stanza ) {
            if ( $line =~ /^\s*s\s+(\d+)/ ) {
                printf {$fh} "  s %d\n", $1;
            }
            elsif ( $line =~ /^\s*b\s+(\d+)\s+(\d+)/ ) {
                printf {$fh} "  b %d %d\n", $S0->($1), $M1->($2);
            }
            elsif ( $line =~ /^\s*e\s+(\d+)\s+(\d+)/ ) {
                printf {$fh} "  e %d %d\n", $S0->($1), $M1->($2);
            }
            elsif ( $line =~ /^\s*l\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ )
            {
                printf {$fh} "  l %d %d %d %d %g\n", $S0->($1), $M1->($2),
                    $S0->($3), $M1->($4), $5;
            }
        }
        print {$fh} "}\n";
    }

    #----------------------------#
    # x-stanza
    #----------------------------#
    if ( $lav =~ /x {\n(.+?)}/s ) {
        my $x_stanza = $1;
        print {$fh} "x {\n$x_stanza}\n";
    }

    #----------------------------#
    # m-stanza
    #----------------------------#
    if ( $lav =~ /m {\n(.+?)}/s ) {
        my $m_stanza = $1;
        print {$fh} "m {\n";
        for my $line ( split /\n/, $m_stanza ) {
            if ( $line =~ /^\s*n\s+(\d+)/ ) {
                printf {$fh} "  n %d\n", $1;
            }
            elsif ( $line =~ /^\s*x\s+(\d+)\s+(\d+)/ ) {
                printf {$fh} "  x %d %d\n", $1 + $t_from - 1,
                    $2 + $t_from - 1;
            }
        }
        print {$fh} "}\n";
    }
}

print {$fh} "#:eof\n";
close $fh;

exit;

__END__

