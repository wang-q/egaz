#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use IO::Zlib;
use App::Fasops::Common qw(decode_header revcom);
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

check_header.pl - check genome location in (blocked) fasta headers

=head1 SYNOPSIS

    perl check_header.pl -i <fasta file> -g <genome file>
      Options:
        --help          -?          brief help message
        --input         -i  STR     be checked file, normal multi fasta or blocked fasta
        --genome        -g  STR     one multi fasta file contains genome sequences
        --name          -n  STR     which species to be checked, omit this will check all sequences
        --detail                    write a fasta file report error sequences

find ~/data/alignment/yeast_genome/S288c/ -name "*.fasta" \
    | sort \
    | xargs cat > ~/data/alignment/yeast_genome/S288c.fasta

perl check_header.pl --in ~/data/alignment/self/yeast_new/S288cvsselfalign_fasta/I.net.axt.gz.fas \
    -g ~/data/alignment/yeast_genome/S288c.fasta \
    --detail

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'input|i=s'  => \my $in_file,
    'genome|g=s' => \my $genome,
    'name|n=s'   => \my $name,
    'detail'     => \my $detail,
) or HelpMessage(1);

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Check headers for [$in_file]");

my $in_fh = IO::Zlib->new( $in_file, "rb" );

my $log_file;
if ($detail) {
    $log_file = path($in_file)->basename . '.log.txt';
    path($log_file)->remove;
}

{
    my $header;
    my $content = '';
    while ( my $line = <$in_fh> ) {
        chomp $line;

        if ( $line =~ /^\>[\w:-]+/ ) {

            # the first sequence is ready
            if ( defined $header ) {
                check_seq( $header, $content, $log_file );
            }

            # prepare to accept next sequence
            $line =~ s/^\>//;
            $header = $line;

            # clean previous sequence
            $content = '';
        }
        elsif ( $line =~ /^[\w-]+/ ) {
            $line =~ s/[^\w]//g;    # Delete '-'s
            $line = uc $line;
            $content .= $line;
        }
        else {                      # Blank line, do nothing
        }
    }

    # for last sequece
    check_seq( $header, $content, $log_file );
}

$in_fh->close;

$stopwatch->end_message( "All sequences scanned.", "duration" );

exit;

sub check_seq {
    my $header   = shift;
    my $seq      = shift;
    my $log_file = shift;

    my $info = decode_header($header);
    if ( $name and $name ne $info->{name} ) {
        return;
    }

    if ( $info->{chr_strand} eq '-' ) {
        $seq = revcom($seq);
    }

    my $location;
    if ( $info->{chr_end} ) {
        $location = sprintf "%s:%s-%s", $info->{chr_name}, $info->{chr_start}, $info->{chr_end};
    }
    else {
        $location = sprintf "%s:%s", $info->{chr_name}, $info->{chr_start};
    }
    my $seq_genome = uc get_seq_faidx( $genome, $location );

    if ( $seq ne $seq_genome ) {
        printf "FAILED\t%s\n", $header;
        if ($log_file) {
            my $str = ">$header\n";
            $str .= "$seq\n";
            $str .= ">$location\n";
            $str .= "$seq_genome\n";
            $str .= "\n";
            path($log_file)->append($str);
        }
    }

    return;
}

sub get_seq_faidx {
    my $genome   = shift;
    my $location = shift;    # I:1-100

    my $cmd = sprintf "samtools faidx %s %s", $genome, $location;
    open my $fh_pipe, '-|', $cmd;

    my $seq;
    while ( my $line = <$fh_pipe> ) {
        chomp $line;
        if ( $line =~ /^[\w-]+/ ) {
            $seq .= $line;
        }
    }
    close($fh_pipe);

    return $seq;
}

__END__
