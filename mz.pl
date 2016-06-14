#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use MCE;
use MCE::Flow;

use Bio::Phylo::IO qw(parse);
use Set::Scalar;
use List::MoreUtils qw(uniq zip);
use Time::Duration;
use Path::Tiny;
use String::Compare;
use Number::Format qw(format_bytes);

use File::Find::Rule;
use File::Copy::Recursive qw(fcopy);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

mz.pl - Multiz step by step

=head1 SYNOPSIS

    perl mz.pl -d <axt dir> -d <axt dir> --tree <.nwk> [options]
      Options:
        --help          -?          brief help message
        --dir           -d  @STR    dirs of *pairwise* maf
        --tree              STR     newick tree, must be a rooted tree
        --target            STR     target name
        --out           -o  STR     output dir
        --all                       don't drop unused syntenies
        --noclean                   keep intermediate file
        --syn                       .synNet.maf or .net.maf
        --parallel      -p  INT     run in parallel mode, [1]


        perl mz.pl  -d ~/data/alignment/arabidopsis19/AthvsLyrata/ \
                    -d ~/data/alignment/arabidopsis19/AthvsBur_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsZu_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsNo_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsLer_0/ \
                    --tree ~/data/alignment/arabidopsis19/19way.nwk \
                    --parallel 8 \
                    -syn

=cut

GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'dir|d=s'  => \my @dirs,
    'out|o=s'  => \my $out_dir,
    'tree=s'   => \my $tree_file,
    'target=s' => \my $target_name,
    'all'      => \my $all,
    'noclean'  => \my $noclean,
    'parallel|p=i' => \( my $parallel = 1 ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

my $suffix = '.maf';
my @files  = File::Find::Rule->file->name("*$suffix")->in(@dirs);
printf "\n----%4s $suffix files ----\n", scalar @files;
if ( scalar @files == 0 ) {
    @files = sort File::Find::Rule->file->name("*$suffix.gz")->in(@dirs);
    printf "\n----Total $suffix.gz Files: %4s----\n\n", scalar @files;
    $suffix = '.maf.gz';
}

if ( scalar @files == 0 ) {
    die "Can't find .maf or .maf.gz files\n";
}

#----------------------------------------------------------#
# Gather species list from maf files
#----------------------------------------------------------#
#---
#Bur_0:
#  1: ~/data/alignment/arabidopsis19/AthvsBur_0/mafNet/chr1.net.maf
#  2: ...
#lyrata_65:
#  1: ~/data/alignment/arabidopsis19/AthvsLyrata_set01_4/mafNet/chr1.net.maf
#  2: ...

my $file_of = {};    # all info here
my %seen;            # count
my @potential_targets;
my @species;         # species list gathered from maf files; and then shift target out

{
    print "Get species list\n";
    my $worker = sub {
        my ( $self, $chunk_ref, $chunk_id ) = @_;
        my $file = $chunk_ref->[0];

        my $cmd
            = "gzip -dcf $file " . q{| perl -nl -e '/^s (\w+)/ or next; print $1' | sort | uniq};
        my @list = grep { defined $_ } split /\n/, `$cmd`;
        if ( @list > 2 ) {
            print "There are three or more species in [$file].\n";
            print YAML::Syck::Dump \@list;
            die;
        }

        MCE->gather( $file, [@list] );
    };
    MCE::Flow::init {
        chunk_size  => 1,
        max_workers => $parallel,
    };
    my %list_of = mce_flow $worker, \@files;
    MCE::Flow::finish;

    print "Assign files to species\n";
    for my $file (@files) {
        my @list = @{ $list_of{$file} };
        $seen{$_}++ for @list;

        my $chr_name = path($file)->basename( ".net$suffix", ".synNet$suffix", $suffix );
        for my $sp (@list) {
            if ( exists $file_of->{$sp} ) {
                if ( exists $file_of->{$sp}{$chr_name} ) {
                    push @potential_targets, $sp;
                    push @{ $file_of->{$sp}{$chr_name} }, $file;
                }
                else {
                    $file_of->{$sp}{$chr_name} = [$file];
                }
                $file_of->{$sp}{chr_set}->insert($chr_name);
            }
            else {
                my $chr_set = Set::Scalar->new;
                $chr_set->insert($chr_name);
                $file_of->{$sp} = { $chr_name => [$file], chr_set => $chr_set };
            }
        }
    }

    @potential_targets = uniq(@potential_targets);
    (@species) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, $seen{$_} ] }
        keys %seen;

    print "\n";
}

#----------------------------#
# check number of species
#----------------------------#
if ( scalar @species < 3 ) {
    print "There're too few species, [@species].\nCan't run multiz.\n";
    print "Just copy pairwise maf to [$out_dir]\n";

    for my $file (@files) {
        fcopy( $file, $out_dir );
    }

    print "\n";
    exit;
}

#----------------------------#
# check target
#----------------------------#
if ( @potential_targets > 1 ) {
    print "There are more than 1 potential targets\n";
    print YAML::Syck::Dump { potential_targets => \@potential_targets };
    die;
}
else {
    if ( !defined $target_name ) {
        $target_name = shift @species;
        if ( $target_name eq $potential_targets[0] ) {
            printf "%s appears %d times, use it as target.\n", $target_name, $seen{$target_name};
        }
        else {
            print "Can't find target.\n";
            die;
        }
    }
    else {
        my ($dummy) = grep { $_ eq $target_name } @species;
        if ( !defined $dummy ) {
            die "Can't find target_name [$target_name] in .maf files.\n";
        }
        elsif ( $target_name eq $potential_targets[0] ) {
            print "Your assigned target [$target_name] is OK\n";
            @species = grep { $_ ne $target_name } @species;
        }
        else {
            print "Your assigned target [$target_name] isn't OK.\n";
            print "It should be [$potential_targets[0]].\n";
            die;
        }
    }
    print "\n";
}

#----------------------------#
# Find chromosomes to be processed
#----------------------------#
my @chr_names = sort $file_of->{$target_name}{chr_set}->members;    # all target chromosomes
{
    print "Target chromosomes are [@chr_names]\n";

    # check other species occurrence number
    my @occurrence = sort { $b <=> $a } uniq( @seen{@species} );
    if ( @occurrence > 1 ) {
        print "Species occurrence number inconsistency [@occurrence]\n";
        print "We will skip some chromosomes\n";
        print YAML::Syck::Dump \%seen;
        print "\n";

        my $intersect_chr_set = $file_of->{$target_name}{chr_set}->clone;
        for my $sp ( keys %{$file_of} ) {
            $intersect_chr_set = $intersect_chr_set->intersection( $file_of->{$sp}{chr_set} );
        }
        @chr_names = sort $intersect_chr_set->members;
        print "Chromosomes to be processed are [@chr_names]\n";
    }
    print "\n";

    # sort @species by distances in tree
    my $ladder = ladder( $tree_file, $target_name );
    my @order = map { ref eq 'ARRAY' ? @$_ : $_ } @{$ladder};
    my %item_order = map { $order[$_] => $_ } 0 .. $#order;

    @species = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { [ $_, $item_order{$_} ] } @species;

    print "Order of stitch [@species]\n";
    print "\n";
}

unless ($out_dir) {
    $out_dir = ucfirst $target_name . "_n" . ( scalars(@species) + 1 );
}
path($out_dir)->mkpath;

YAML::Syck::DumpFile(
    path( $out_dir, 'info.yml' )->stringify,
    {   file_of   => $file_of,
        chr_names => \@chr_names,
    }
);

#----------------------------------------------------------#
# Finally, multiz comes
#----------------------------------------------------------#

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;

    my $chr_name = $chunk_ref->[0];

# multiz.v11.2: -- aligning two files of alignment blocks where top rows are
# always the reference, reference in both files cannot have duplicats
# args: [R=?] [M=?] file1 file2 v? [out1 out2] [nohead] [all]
#         R(30) radius in dynamic programming.
#         M(1) minimum output width.
#         out1 out2(null) null: stdout; out1 out2: file names for collecting unused input.
#         nohead(null) null: output maf header; nohead: not to output maf header.
#         all(null) null: not to output single-row blocks; all: output all blocks.
#
# multiz ~/data/alignment/arabidopsis19/AthvsZu_0/mafSynNet/chr1.synNet.maf \
#       ~/data/alignment/arabidopsis19/AthvsNo_0/mafSynNet/chr1.synNet.maf \
#       1 out1 out2 > step1.chr1.maf
# multiz step1.chr1.maf ~/data/alignment/arabidopsis19/AthvsBur_0/mafSynNet/chr1.synNet.maf \
#       1 out1 out2 > step2.chr1.maf
# multiz step2.chr1.maf ~/data/alignment/arabidopsis19/AthvsLyrata_set01_4/mafSynNet/chr1.synNet.maf
#       1 out1 out2 > step3.chr1.maf

    my @species_copy = @species;
    my ( $species1, $species2 );
    my $maf_step;
    my $step = 1;
    my $str  = '';
    while (@species_copy) {
        my ( $maf1, $maf2 );

        if ( !defined $species1 ) {
            $species1 = shift @species_copy;
            $maf1     = $file_of->{$species1}{$chr_name}[0];
        }
        else {
            $maf1 = $maf_step;
        }
        if ( !defined $species2 ) {
            $species2 = shift @species_copy;
            $maf2     = $file_of->{$species2}{$chr_name}[0];
        }

        my $out1 = "$out_dir/$chr_name.out1";
        my $out2 = "$out_dir/$chr_name.out2";

        $maf_step
            = @species_copy
            ? "$out_dir/$chr_name.step$step.maf"
            : "$out_dir/$chr_name.maf";

        # here we set out1 and out2 to discard unused synteny
        # Omit out1 and out2, unused synteny will be printed to stdout and
        # reused by following multiz processes
        print "Run multiz...\n";
        my $cmd
            = "multiz" . " M=10"
            . " $maf1"
            . " $maf2" . " 1 "
            . " $out1"
            . " $out2"
            . " > $maf_step";
        exec_cmd($cmd);
        print "Step [$step] .maf file generated.\n\n";

        $str .= "$chr_name.step$step,";
        $str .= "$species1,$species2,";
        for my $file ( $maf1, $maf2, $out1, $out2, $maf_step ) {
            $str .= format_bytes( -s $file, base => 1000 );
            $str .= ",";
        }
        $str .= format_bytes( ( -s $maf_step ) / ( $step + 2 ), base => 1000 );
        $str .= "\n";

        $species1 = "step$step";
        $species2 = undef;
        $step++;
    }

    if ( !$noclean ) {
        print "Clean temp files.\n";
        path( $out_dir, "$chr_name.out1" )->remove;
        path( $out_dir, "$chr_name.out2" )->remove;
        for ( path($out_dir)->children(qr/^$chr_name\.step/) ) {
            $_->remove;
        }
    }

    my $cmd = "gzip " . "$out_dir/$chr_name.maf";
    exec_cmd($cmd);

    print $str;
    path( $out_dir, "$chr_name.temp.csv" )->remove;
    path( $out_dir, "$chr_name.temp.csv" )->spew($str);

    return;
};

my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
$mce->foreach( \@chr_names, $worker );

{    # summary
    my $cmd = "echo 'step,spe1,spe2,maf1,maf2,out1,out2,size,per_size'" . " > $out_dir/steps.csv";
    exec_cmd($cmd);

    $cmd
        = "find $out_dir -type f -name '*.temp.csv'"
        . " | sort -n"
        . " | xargs cat >> $out_dir/steps.csv";
    exec_cmd($cmd);

    $cmd = "find $out_dir -type f -name '*.temp.csv'" . " | xargs rm";
    exec_cmd($cmd);
}

print "\n";
print "Runtime ", duration( time - $start_time ), ".\n";
print "=" x 30, "\n";

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

#----------------------------#
# parsing tree
#----------------------------#
sub ladder {
    my $file        = shift;
    my $target_name = shift;

    my $tree = parse( -format => 'newick', -file => $file )->first;

    my $target = $tree->get_by_name($target_name);

    my @ladder = ( [$target_name], );
    my $set = Set::Scalar->new;
    $set->insert($target_name);

    my $start = $target;
    while (1) {
        my $parent = $start->get_parent;
        last unless $parent;
        my @nodes
            = grep { !$set->has( $_->get_name ) } @{ $parent->get_terminals };
        $set->insert( $_->get_name ) for @nodes;
        my $distance_of = {};
        for my $node (@nodes) {
            $distance_of->{ $node->get_name } = $node->calc_patristic_distance($target);
        }
        my @sorted = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { [ $_, $distance_of->{$_} ] }
            keys %{$distance_of};
        push @ladder, [@sorted];

        $start = $parent;
    }

    return \@ladder;
}

sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

__END__
