#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(any all uniq);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $op = 'intersect';

my $file1;
my $file2;

my $outfile;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'op|operation=s' => \$op,
    'f1|file1=s' => \$file1,
    'f2|file2=s' => \$file2,
    'o|outfile=s' => \$outfile,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

if ( !-e $file1 or !-e $file2 ) {
    print "You should provide 2 YAML files\n";
    exit;
}

if ( $op =~ /^dif/i ) {
    $op = 'diff';
}
elsif ( $op =~ /^uni/i ) {
    $op = 'union';
}
elsif ( $op =~ /^int/i ) {
    $op = 'intersect';
}

unless ($outfile) {
    $outfile = "$op.yml";
}

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Compare runlist...");

print "Loading...\n";
my @sets;
my %seen;
for my $file ( $file1, $file2 ) {
    print "Filename: $file\n";
    my $runlist_of = LoadFile($file);

    my $set_of = {};
    for my $key ( sort keys %{$runlist_of} ) {
        $seen{$key}++;
        my $set = AlignDB::IntSpan->new( $runlist_of->{$key} );
        printf "key:\t%s\tlength:\t%s\n", $key, $set->size;
        $set_of->{$key} = $set;
    }
    push @sets, $set_of;
}
print "\n";

my @keys = grep {$seen{$_} >= 2} sort keys %seen;
my $op_runlist_of = {};
for my $key (@keys) {
    print "For $key\n";
    my $op_set = $sets[0]->{$key}->$op($sets[1]->{$key});
    next if $op_set->is_empty;
    $op_runlist_of->{$key} = $op_set->runlist;
    print " " x 4, $op_set->size, "\n";
}

DumpFile( $outfile, $op_runlist_of );

$stopwatch->end_message;

__END__

=head1 NAME

    compare_runlist.pl - compare 2 chromosome runlists

=head1 SYNOPSIS

    perl compare_runlist.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --op                operations: intersect, union, diff
        --file1
        --file2
        --outfile

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

