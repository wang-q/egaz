#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(exec_cmd);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

fasta_blastn.pl - Blasting between two fasta files
    
=head1 SYNOPSIS

    perl fasta_blastn.pl -f <fasta file> -g <genome file> [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     query fasta file
        --genome        -g  STR     reference genome file
        --output        -o  STR     output
        --parallel      -p  INT     default is [8]

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \( my $file ),
    'genome|g=s' => \( my $genome ),
    'output|o=s' => \( my $output ),
    'parallel|p=i' => \( my $parallel = 8 ),
) or HelpMessage(1);

if ( !defined $file ) {
    die "Need --file\n";
}
elsif ( !path($file)->is_file ) {
    die "--file [$file] doesn't exist\n";
}

if ( !defined $genome ) {
    die "--genome is needed\n";
}
elsif ( !path($genome)->is_file ) {
    die "--genome doesn't exist\n";
}

if ( !$output ) {
    $output = path($file)->basename;
    ($output) = grep {defined} split /\./, $output;
    $output = "$output.blast";
}
path($output)->remove;

#----------------------------------------------------------#
# Run
#----------------------------------------------------------#
$stopwatch->start_message("Start...");

$stopwatch->block_message("Build blast db");
exec_cmd("makeblastdb -dbtype nucl -in $genome");

$stopwatch->block_message("Blasting...");
my $cmd
    = sprintf "blastn -task megablast -evalue 0.0001 -word_size 40"    # megablast with word size 40
    . " -max_target_seqs 10 -max_hsps 5 -culling_limit 10"             # reduce size of reports
    . " -dust no -soft_masking false"
    . " -outfmt '7 qseqid sseqid qstart qend sstart send qlen slen nident'"
    . " -num_threads %d -db %s -query %s"
    . " > %s", $parallel, $genome, $file, $output;
exec_cmd($cmd);

$stopwatch->end_message;

exit;

__END__
