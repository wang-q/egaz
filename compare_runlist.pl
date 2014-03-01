#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(any all uniq);
use Set::Scalar;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $op = 'intersect';

my $file1;
my $file2;

my $multi_key;    # file1 has multiple keys

my $outfile;

my $remove_chr;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'         => \$help,
    'man'            => \$man,
    'op|operation=s' => \$op,
    'f1|file1=s'     => \$file1,
    'f2|file2=s'     => \$file2,
    'mk'             => \$multi_key,
    'o|outfile=s'    => \$outfile,
    'r|remove'       => \$remove_chr,
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
elsif ( $op =~ /^xor/i ) {
    $op = 'xor';
}

unless ($outfile) {
    $outfile = "$op.yml";
}

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Compare runlist...");

#----------------------------#
# Loading
#----------------------------#
my $all_name_set = Set::Scalar->new;

# file1
print "Loading $file1...\n";
my $s1_of = {};
my @keys;
if ($multi_key) {
    my $yml = LoadFile($file1);
    @keys = sort keys %{$yml};

    for my $key (@keys) {
        my ( $chr_set_of, $chr_name_set )
            = runlist2set( $yml->{$key}, $remove_chr );
        $s1_of->{$key} = $chr_set_of;
        $all_name_set->insert( $chr_name_set->members );
    }
}
else {
    @keys = ("__single");

    my ( $chr_set_of, $chr_name_set )
        = runlist2set( LoadFile($file1), $remove_chr );
    $s1_of->{__single} = $chr_set_of;
    $all_name_set->insert( $chr_name_set->members );
}

# file2
print "Loading $file2...\n";
my $s2;
{
    my ( $chr_set_of, $chr_name_set )
        = runlist2set( LoadFile($file2), $remove_chr );
    $s2 = $chr_set_of;
    $all_name_set->insert( $chr_name_set->members );
}

#----------------------------#
# Operation
#----------------------------#
print "\nOperation $op\n";
my $op_result_of = { map { $_ => {} } @keys };

for my $key (@keys) {
    print "For key $key\n";
    my $s1 = $s1_of->{$key};

    # give empty set to non-existing chrs
    for my $s ( $s1, $s2 ) {
        for my $chr ( $all_name_set->members ) {
            if ( !exists $s->{$chr} ) {
                $s->{$chr} = AlignDB::IntSpan->new;
            }
        }
    }

    # operate on each chr
    for my $chr ( sort $all_name_set->members ) {
        print "\tFor chr $chr\n";
        my $op_set = $s1->{$chr}->$op( $s2->{$chr} );
        $op_result_of->{$key}{$chr} = $op_set->runlist;
        print "\t\tlength:\t", $op_set->size, "\n";
    }
}

if ($multi_key) {
    DumpFile( $outfile, $op_result_of );
}
else {
    DumpFile( $outfile, $op_result_of->{__single} );
}

$stopwatch->end_message;

sub runlist2set {
    my $runlist_of = shift;
    my $remove_chr = shift;

    my $set_of       = {};
    my $chr_name_set = Set::Scalar->new;

    for my $chr ( sort keys %{$runlist_of} ) {
        my $new_chr = $chr;
        $new_chr =~ s/chr0?// if $remove_chr;
        $chr_name_set->insert($chr);
        my $set = AlignDB::IntSpan->new( $runlist_of->{$chr} );
        printf "\tkey:\t%s\tlength:\t%s\n", $new_chr, $set->size;
        $set_of->{$new_chr} = $set;
    }

    return ( $set_of, $chr_name_set );
}

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

