#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Basename;
use File::Find::Rule;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $dir;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'  => \$help,
    'man'     => \$man,
    'd|dir=s' => \$dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Merge chr runlist...");

if ( !-d $dir ) {
    print "You should provide a dir\n";
    exit;
}

my @yml_files = File::Find::Rule->file->name( '*.yaml', '*.yml' )->in($dir);
printf "\n----Total YAML Files: %4s----\n\n", scalar @yml_files;

my $master = {};
my @dirs = grep {$_} File::Spec->splitdir($dir);
push @dirs, $dirs[-1] . ".yml";
my $file_out = File::Spec->catfile(@dirs);

for my $file (@yml_files) {
    my ( $base, $dir ) = fileparse( $file, ".yaml", ".yml" );
    my ($word) = split /[^\w]+/, $base;

    my $content = LoadFile($file);
    $master->{$word} = $content;

    unlink $file;
}
DumpFile( $file_out, $master );

open my $report_fh, '>', "$file_out.txt";
for my $key (sort keys %{$master}) {
    my $set = AlignDB::IntSpan->new($master->{$key});
    printf {$report_fh} "key:\t[%s]\tlength:\t[%s]\n", $key, $set->size;
}
close $report_fh;

$stopwatch->end_message;

__END__

=head1 NAME

    write_axt.pl - extract sequence of a certain feature from alignDB

=head1 SYNOPSIS

    write_axt.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --ensembl           ensembl database name

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

