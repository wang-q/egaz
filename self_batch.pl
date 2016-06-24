#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use IPC::Cmd qw(can_run);
use Path::Tiny;
use Template;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
);

=head1 NAME

self_batch.pl - Full procedure for self genome alignments.

=head1 SYNOPSIS

    perl self_batch.pl [options]
      Options:
        --help          -?          brief help message
        --working_dir   -w  STR     Default is [.]
        --seq_dir       -s  STR     Will do prep_fa() from this dir or use seqs store in $working_dir
        --target        -t  STR
        --queries       -q  @STR
        --csv_taxon     -c  STR     All taxons in this project (may also contain unused taxons)
        --length            INT     Minimal length of paralogous fragments
        --name_str      -n  STR     Default is []
        --parted                    Sequences are partitioned
        --noblast                   Don't blast against genomes
        --msa               STR     Aligning program for refine. Default is [mafft]
        --norm                      RepeatMasker has been done.
        --nostat                    Don't do stat stuffs
        --norawphylo                Skip rawphylo
        --parallel          INT     number of child processes

=cut

my $aligndb = path( $FindBin::RealBin, "..", "alignDB" )->absolute->stringify;
my $egaz = path($FindBin::RealBin)->absolute->stringify;

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'working_dir|w=s' => \( my $working_dir = "." ),
    'seq_dir|s=s'     => \my $seq_dir,
    'target|t=s'      => \my $target,
    'queries|q=s'     => \my @queries,
    'csv_taxon|c=s'   => \my $csv_taxon,
    'length=i'        => \( my $length      = 1000 ),
    'name_str|n=s'    => \( my $name_str    = "working" ),
    'parted'          => \my $parted,
    'noblast'         => \my $noblast,
    'msa=s'           => \( my $msa         = 'mafft' ),
    'norm'            => \my $norm,
    'nostat'          => \my $nostat,
    'parallel=i'      => \( my $parallel    = 4 ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains summary...");

# prepare working dir
{
    print "Working on $name_str\n";
    $working_dir = path( $working_dir, $name_str )->absolute;
    $working_dir->mkpath;
    $working_dir = $working_dir->stringify;
    print " " x 4, "Working dir is $working_dir\n";

    path( $working_dir, 'Genomes' )->mkpath;
    path( $working_dir, 'Pairwise' )->mkpath;
    path( $working_dir, 'Processing' )->mkpath;
    path( $working_dir, 'Results' )->mkpath;
}

# build basic information
my @data;
{
    my %taxon_of;
    if ($csv_taxon) {
        for my $line ( path($csv_taxon)->lines ) {
            my @fields = split /,/, $line;
            if ( $#fields >= 2 ) {
                $taxon_of{ $fields[0] } = $fields[1];
            }
        }
    }
    @data = map {
        {   name  => $_,
            taxon => exists $taxon_of{$_} ? $taxon_of{$_} : 0,
            dir   => path( $working_dir, 'Genomes', $_ )->stringify,
        }
    } ( $target, @queries );
}

# if seqs is not in working dir, copy them from seq_dir
if ($seq_dir) {
    print "Get seqs from [$seq_dir]\n";

    for my $id ( $target, @queries ) {
        print " " x 4 . "Copy seq of [$id]\n";

        my $original_dir = path( $seq_dir, $id )->stringify;
        my $cur_dir = path( $working_dir, 'Genomes', $id );
        $cur_dir->mkpath;
        $cur_dir = $cur_dir->stringify;

        my @fa_files
            = File::Find::Rule->file->name( '*.fna', '*.fa', '*.fas',
            '*.fasta' )->in($original_dir);

        printf " " x 8 . "Total %d fasta file(s)\n", scalar @fa_files;

        for my $fa_file (@fa_files) {
            my $basename = prep_fa( $fa_file, $cur_dir );

            my $gff_file = path( $original_dir, "$basename.gff" );
            if ( $gff_file->is_file ) {
                $gff_file->copy($cur_dir);
            }
            my $rm_gff_file = path( $original_dir, "$basename.rm.gff" );
            if ( $rm_gff_file->is_file ) {
                $rm_gff_file->copy($cur_dir);
            }
        }
    }
}

{
    my $tt = Template->new( ABSOLUTE => 1, );
    my $sh_name;

    #----------------------------#
    # all *.sh files
    #----------------------------#

    # real_chr.sh
    $sh_name = "1_real_chr.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "1_real_chr.tt2" )->stringify,
        {   data        => \@data,
            working_dir => $working_dir,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # file-rm.sh
    if ( !$norm ) {
        $sh_name = "2_file_rm.sh";
        print "Create $sh_name\n";
        $tt->process(
            path( $FindBin::RealBin, "template", "2_file_rm.tt2" )->stringify,
            {   data        => \@data,
                parallel    => $parallel,
                working_dir => $working_dir,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # self_cmd.sh
    $sh_name = "3_self_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "3_self_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            parted      => $parted,
            name_str    => $name_str,
            all_ids     => [ $target, @queries ],
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # proc_cmd.sh
    $sh_name = "4_proc_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "4_proc_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            msa         => $msa,
            noblast     => $noblast,
            name_str    => $name_str,
            all_ids     => [ $target, @queries ],
            data        => \@data,
            length      => $length,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # circos.conf
    for my $id ( $target, @queries ) {
        print "    Create circos.conf for $id\n";
        $tt->process(
            path( $FindBin::RealBin, "template", "circos.tt2" )->stringify,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                name_str    => $name_str,
                id          => $id,
            },
            path( $working_dir, 'Processing', "${id}", "circos.conf" )
                ->stringify
        ) or die Template->error;
    }

    # circos_cmd.sh
    $sh_name = "5_circos_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "5_circos_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            name_str    => $name_str,
            all_ids     => [ $target, @queries ],
            data        => \@data,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # feature_cmd.sh
    $sh_name = "6_feature_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "6_feature_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            name_str    => $name_str,
            all_ids     => [ $target, @queries ],
            data        => \@data,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    if ( !$nostat ) {
        $sh_name = "7_pair_stat.sh";
        print "Create $sh_name\n";
        $tt->process(
            path( $FindBin::RealBin, "template", "7_pair_stat.tt2" )->stringify,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                name_str    => $name_str,
                all_ids     => [ $target, @queries ],
                data        => \@data,
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # pack_it_up.sh
    $sh_name = "9_pack_it_up.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "9_pack_it_up.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            name_str    => $name_str,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # message
    $stopwatch->block_message("Execute *.sh files in order.");
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub prep_fa {
    my $infile = shift;
    my $dir    = shift;

    my $basename = path($infile)->basename( '.fna', '.fa', '.fas', '.fasta' );
    my $in_fh    = path($infile)->openr;
    my $out_fh   = path( $dir, "$basename.fa" )->openw;
    while (<$in_fh>) {
        if (/>/) {
            print {$out_fh} ">$basename\n";
        }
        else {
            print {$out_fh} $_;
        }
    }
    close $out_fh;
    close $in_fh;

    return $basename;
}

__END__
