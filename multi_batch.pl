#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use DBI;
use Text::CSV_XS;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);
use Template;

use IPC::Cmd qw(can_run);
use Path::Tiny;
use File::Find::Rule;

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

multi_batch.pl - Full procedure for multiple genome alignments.

=head1 SYNOPSIS

    perl multi_batch.pl [options]
      Options:
        --help          -?          brief help message
        --working_dir   -w  STR     Default is [.]
        --seq_dir       -s  STR     Will do prep_fa() from this dir or use seqs store in $working_dir
        --target        -t  STR
        --queries       -q  @STR
        --outgroup      -o  STR
        --csv_taxon     -c  STR     All taxons in this project (may also contain unused taxons)
        --length            INT     Minimal length of orthologous fragments
        --name_str      -n  STR     Default is [working]
        --phylo_tree    -p  STR     Predefined phylogenetic tree
        --multi_name    -m  STR     Naming multiply alignment, the default value is $name_str
                                    This option is for more than one align combination.
        --msa               STR     Aligning program for refine. Default is [mafft]
        --norm                      RepeatMasker has been done.
        --nostat                    Don't do stat stuffs
        --norawphylo                Skip rawphylo
        --parallel          INT     Number of child processes. Default is [4]

=cut

my $aligndb = path( $FindBin::RealBin, "..", "alignDB" )->absolute->stringify;
my $egaz = path($FindBin::RealBin)->absolute->stringify;

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'working_dir|w=s' => \( my $working_dir = "." ),
    'seq_dir|s=s'     => \my $seq_dir,
    'target|t=s'      => \my $target,
    'queries|q=s'     => \my @queries,
    'outgroup|o|r=s'  => \my $outgroup,
    'csv_taxon|c=s'   => \my $csv_taxon,
    'length=i'        => \( my $length      = 1000 ),
    'name_str|n=s'    => \( my $name_str    = "working" ),
    'phylo_tree|p=s'  => \my $phylo_tree,
    'multi_name|m=s'  => \my $multi_name,
    'msa=s'           => \( my $msa         = 'mafft' ),
    'norm'            => \my $norm,
    'nostat'          => \my $nostat,
    'norawphylo'      => \my $norawphylo,
    'parallel=i'      => \( my $parallel    = 4 ),
) or Getopt::Long::HelpMessage(1);

if ( defined $phylo_tree ) {
    if ( !-e $phylo_tree ) {
        warn "$phylo_tree does not exists. Unset it.\n";
        undef $phylo_tree;
    }
}

if ( !defined $multi_name ) {
    $multi_name = $name_str;
}

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
    path( $working_dir, 'Stats' )->mkpath;
}

# move $outgroup to last
if ($outgroup) {
    my ($exist) = grep { $_ eq $outgroup } @queries;
    if ( !defined $exist ) {
        die "outgroup does not exist!\n";
    }

    @queries = grep { $_ ne $outgroup } @queries;
    push @queries, $outgroup;
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

# If there's no phylo tree, generate a fake one.
if ( !defined $phylo_tree ) {
    print "Create fake_tree.nwk\n";
    my $fh = path( $working_dir, "fake_tree.nwk" )->openw;
    print {$fh} "(" x scalar(@queries) . "$target";
    for my $id (@queries) {
        print {$fh} ",$id)";
    }
    print {$fh} ";\n";
    close $fh;
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

    # pair_cmd.sh
    $sh_name = "3_pair_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "3_pair_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target      => $target,
            queries     => \@queries,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    # rawphylo.sh
    if ( !$norawphylo and !defined $phylo_tree ) {
        $sh_name = "4_rawphylo.sh";
        print "Create $sh_name\n";
        $tt->process(
            path( $FindBin::RealBin, "template", "4_rawphylo.tt2" )->stringify,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                egaz        => $egaz,
                target      => $target,
                outgroup    => $outgroup,
                queries     => \@queries,
                length      => $length,
                multi_name  => $multi_name,
                avx         => can_run('raxmlHPC-PTHREADS-AVX'),
            },
            path( $working_dir, $sh_name )->stringify
        ) or die Template->error;
    }

    # multi_cmd.sh
    $sh_name = "5_multi_cmd.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "5_multi_cmd.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target      => $target,
            outgroup    => $outgroup,
            queries     => \@queries,
            phylo_tree  => $phylo_tree,
            multi_name  => $multi_name,
            msa         => $msa,
            avx         => can_run('raxmlHPC-PTHREADS-AVX'),
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

# var_list.sh
# Fixme: The column names do not match; the column "Ace_pasteurianus_386B.NC_021991(+):63873-133498" no present in [/tmp/fas_vVWwMHt4/Ace_pasteurianus_IFO_3283_01.NC_013209.+.2789416-2790912.fas.vcf].
    $sh_name = "6_var_list.sh";
    print "Create $sh_name\n";
    $tt->process(
        path( $FindBin::RealBin, "template", "6_var_list.tt2" )->stringify,
        {   stopwatch   => $stopwatch,
            parallel    => $parallel,
            working_dir => $working_dir,
            egaz        => $egaz,
            target      => $target,
            multi_name  => $multi_name,
        },
        path( $working_dir, $sh_name )->stringify
    ) or die Template->error;

    if ( !$nostat ) {
        $sh_name = "7_multi_db_only.sh";
        print "Create $sh_name\n";
        $tt->process(
            path( $FindBin::RealBin, "template", "7_multi_db_only.tt2" )->stringify,
            {   stopwatch   => $stopwatch,
                parallel    => $parallel,
                working_dir => $working_dir,
                aligndb     => $aligndb,
                target      => $target,
                length      => $length,
                multi_name  => $multi_name,
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
            name_str    => $multi_name,
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
