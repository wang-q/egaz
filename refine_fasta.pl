#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use File::Find::Rule;
use MCE;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

refine_fasta.pl - realign fasta file

=head1 SYNOPSIS

    perl refine_fasta.pl [options]
      Options:
        --help          -?          brief help message
        --input         -i  STR     fasta files' location
        --output        -o  STR     output location
        --msa                       alignment program ([none] means don't do realigning)
        --block                     input is blocked fasta
        --quick                     use quick mode
        --outgroup                  has outgroup at the end
        --expand                    in quick mode, expand indel region
        --join                      in quick mode, join adjacent indel regions
        --parallel                  run in parallel mode

    perl refine_fasta.pl -i G:/S288CvsRM11 --msa muscle --quick

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'input|i=s' => \( my $in_dir = '.' ),
    'output|o=s' => \( my $out_dir ),
    'msa=s'      => \( my $aln_prog = 'clustalw' ),
    'block'      => \my $block,
    'quick'      => \my $quick_mode,
    'outgroup'   => \my $outgroup,
    'expand=i'   => \( my $indel_expand = 50 ),
    'join=i'     => \( my $indel_join = 50 ),
    'parallel=i' => \( my $parallel = 1 ),
) or HelpMessage(1);

#----------------------------------------------------------#
# make output dir
#----------------------------------------------------------#
unless ($out_dir) {
    $out_dir = path($in_dir)->absolute->stringify . "_$aln_prog";
    $out_dir = $out_dir . "_quick" if $quick_mode;
}

if ( path($out_dir)->exists ) {
    warn "$out_dir exists, remove it.\n";
    path($out_dir)->remove_tree;
}
path($out_dir)->mkpath;

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my @files = File::Find::Rule->file->name( '*.fa', '*.fas', '*.fasta' )->in($in_dir);
printf "\n----Total .fas Files: %4s----\n\n", scalar @files;

#----------------------------------------------------------#
# realign
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $infile = $chunk_ref->[0];

    my $stopwatch = AlignDB::Stopwatch->new;
    print "Process $infile\n";

    my ( $seq_of, $seq_names ) = read_fasta($infile);

    if ( $aln_prog ne 'none' ) {
        if ($quick_mode) {
            realign_quick(
                $seq_of,
                $seq_names,
                {   indel_expand => $indel_expand,
                    indel_join   => $indel_join,
                    aln_prog     => $aln_prog,
                }
            );
        }
        else {
            realign_all( $seq_of, $seq_names );
        }
    }

    trim_pure_dash( $seq_of, $seq_names );

    if ($outgroup) {
        trim_outgroup( $seq_of, $seq_names );
    }
    if ($outgroup) {
        trim_complex_indel( $seq_of, $seq_names );
    }

    my $outfile = path($infile)->basename;
    $outfile = $out_dir . "/$outfile";

    open my $out_fh, '>', $outfile;
    for my $name ( @{$seq_names} ) {
        my $seq = $seq_of->{$name};
        print {$out_fh} ">", $name, "\n";
        print {$out_fh} $seq, "\n";
    }
    close $out_fh;
    print "Done.\n\n";
};

my $worker_block = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $infile = $chunk_ref->[0];

    my $stopwatch = AlignDB::Stopwatch->new;
    print "Process $infile\n";

    # don't use $/ = "\n\n", which cause bioperl panic
    open my $in_fh, "<", $infile;
    my $content = '';
    my $count   = 0;
    while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';
            die "headers not equal to seqs\n" if @lines % 2;

            $count++;
            printf "Block [%s] ", $count;
            my ( $seq_of, $seq_names ) = ( {}, [] );

            # store simplified names
            # because names containing . + () |
            my $names = [ 0 .. @lines / 2 - 1 ];
            while (@lines) {
                my $name = shift @lines;
                $name =~ s/^\>//;
                chomp $name;
                my $seq = shift @lines;
                chomp $seq;
                push @{$seq_names}, $name;
                my $idx = scalar @{$seq_names} - 1;
                $seq_of->{$idx} = $seq;
            }

            printf "length [%s]\n", length $seq_of->{0};
            if ( $aln_prog ne 'none' ) {
                if ($quick_mode) {
                    realign_quick(
                        $seq_of, $names,
                        {   indel_expand => $indel_expand,
                            indel_join   => $indel_join,
                            aln_prog     => $aln_prog,
                        }
                    );
                }
                else {
                    realign_all( $seq_of, $names );
                }
            }

            trim_pure_dash( $seq_of, $names );

            if ($outgroup) {
                trim_outgroup( $seq_of, $names );
            }
            if ($outgroup) {
                trim_complex_indel( $seq_of, $names );
            }

            my $outfile = path($infile)->basename;
            $outfile = $out_dir . "/$outfile";

            open my $out_fh, '>>', $outfile;
            for my $i ( @{$names} ) {
                print {$out_fh} ">", $seq_names->[$i], "\n";
                print {$out_fh} $seq_of->{$i}, "\n";
            }
            print {$out_fh} "\n";
            close $out_fh;
        }
        else {
            $content .= $line;
        }
    }
    close $in_fh;
    print "Done.\n\n";
};

# process each .fasta files
my $stopwatch = AlignDB::Stopwatch->new;

my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
$mce->foreach( [ sort @files ], $block ? $worker_block : $worker );

$stopwatch->block_message( "All files have been processed.", "duration" );
exit;

#----------------------------#
# realign all seqs
#----------------------------#
sub realign_all {
    my $seq_of    = shift;
    my $seq_names = shift;

    my @seqs;
    for ( @{$seq_names} ) {
        push @seqs, $seq_of->{$_};
    }

    my $realigned_seqs = multi_align( \@seqs, $aln_prog );

    for my $i ( 0 .. scalar @{$seq_names} - 1 ) {
        $seq_of->{ $seq_names->[$i] } = uc $realigned_seqs->[$i];
    }

    return;
}

__END__
