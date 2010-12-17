#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use WWW::Mechanize;

my $start_page     = "http://hgdownload.cse.ucsc.edu/downloads.html";
my $url_text_regex = "pairwise alignments";
my $content_regex  = "scoring matrix";
my $wanted_file    = "README.txt";
my $out_dir        = "wanted";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'out_dir=s' => \$out_dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

$|++;

# make lav dir
unless (-e $out_dir) {
    mkdir $out_dir, 0777
      or die "Cannot create \"$out_dir\" directory: $!";
}

my $mech = WWW::Mechanize->new();
$mech->proxy([ 'http', 'ftp' ], 'http://127.0.0.1:8088/');

$mech->get($start_page);

my @links = $mech->find_all_links(text_regex => qr/$url_text_regex/);

my ($count_all, $count_keeped) = (0, 0);
foreach my $link (@links) {
    my $url = $link->url_abs;
    $url =~ /(\w+)\/(\w+)\/?$/;
    my $outfile  = $out_dir . "/" . ucfirst $1 . $2 . ".txt";
    my $filename = $url . "/$wanted_file";

    print "Fetching $filename\n";
    $mech->get($filename);

    my $content = $mech->content(format => "text");

    if ($content =~ /$content_regex/) {
        $content =~ s/^.+miller_lab//s;
        $content =~ s/\-{5,}.+$//s;
        
        print "Saving $outfile\n";
        open my $out_fh, ">$outfile";
        print $out_fh $content;
        close $out_fh;
        $count_keeped++;
    }
    
    $count_all++;
    print "\n";
}

print "all files: $count_all\n" . "keeped files:$count_keeped\n";

exit(0);
