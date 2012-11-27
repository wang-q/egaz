#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Run;
use AlignDB::Util qw(:all);

use Bio::Phylo::IO qw(parse);
use Set::Scalar;
use List::Flatten;
use List::MoreUtils qw(uniq zip);
use Time::Duration;
use Roman;
use File::Find::Rule;
use File::Spec;
use File::Basename;
use File::Remove qw(remove);
use Path::Class;
use String::Compare;
use Number::Format qw(format_bytes);

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# executable file location
my $kent_bin   = "~/bin/x86_64";
my $multiz_bin = "~/bin";

# dirs of *pairwise* maf
my @dirs;

# newick tree
# Must be a rooted tree
my $tree_file;

#  names
my $target_name;

# output dir
my $out_dir;

# .synNet.maf or .net.maf
my $syn;

# don't drop unused syntenies
my $all;

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'         => \$help,
    'man'            => \$man,
    'bin|kent_bin=s' => \$kent_bin,
    'multiz_bin=s'   => \$multiz_bin,
    'd|dir=s'        => \@dirs,
    'out=s'          => \$out_dir,
    'tree=s'         => \$tree_file,
    'target=s'       => \$target_name,
    'syn'            => \$syn,
    'all'            => \$all,
    'p|parallel=i'   => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $start_time = time;
print "\n", "=" x 30, "\n";
print "Processing...\n";

my $suffix = $syn ? '.synNet.maf' : '.net.maf';
my $gzip;
my @files = File::Find::Rule->file->name("*$suffix")->in(@dirs);
printf "\n----%4s $suffix files ----\n", scalar @files;
if ( scalar @files == 0 ) {
    @files = sort File::Find::Rule->file->name("*$suffix.gz")->in(@dirs);
    printf "\n----Total $suffix.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

if ( scalar @files == 0 ) {
    print "Can't find maf files\n";
    exit;
}

#----------------------------#
# decompress
#----------------------------#
if ($gzip) {
    print "maf files is gzipped, multiz can't handle them.\n";
    print "Decompress...\n";
    
    my $worker = sub {
        my $job = shift;

        my $file = $job;
        my $cmd  = "gzip -d $file";
        exec_cmd($cmd);

        return;
    };

    my @jobs = sort @files;
    my $run  = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;
    
    # refind newly decompressed maf files
    @files = File::Find::Rule->file->name("*$suffix")->in(@dirs);
}

#----------------------------#
# Gather species list from maf and tree
#----------------------------#
my @species;
my %species_of;
my $number_of_chr;
{
    print "Get species list\n";
    for (@files) {
        my $cmd = "$kent_bin/mafSpeciesList $_ stdout";
        $species_of{$_} = [ grep {$_} split /\n/, `$cmd` ];
    }

    # count
    my %seen;
    for my $key ( keys %species_of ) {
        for my $elem ( @{ $species_of{$key} } ) {
            $seen{$elem}++;
        }
    }
    (@species) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, $seen{$_} ] }
        keys %seen;

    # check species number
    my $number = scalar @species;
    die "There're only $number species, [@species].\nCan't run multiz.\n"
        if $number < 3;

    if ( !$target_name ) {
        printf "%s appears %d times, use it as target.\n", $species[0],
            $seen{ $species[0] };
        $target_name = $species[0];
        shift @species;
    }
    else {
        my ($dummy) = grep { $_ eq $target_name } @species;
        die "Can't find target_name [$target_name] in .maf files.\n"
            unless $dummy;
        @species = grep { $_ ne $target_name } @species;
    }

    # check other species occurrence number
    my @occurrence = uniq( @seen{@species} );
    if ( @occurrence > 1 ) {
        print Dump \%seen;
        die "Species occurrence number inconsistency [@occurrence]\n";
    }

    $number_of_chr = $occurrence[0];
    print "I guess there are $number_of_chr chromosomes.\n";

    # sort @species by distances in tree
    my $ladder     = ladder( $tree_file, $target_name );
    my @order      = flat( @{$ladder} );
    my %item_order = map { $order[$_] => $_ } 0 .. $#order;

    @species = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { [ $_, $item_order{$_} ] } @species;

    print "Order of stitch [@species]\n";

    my $cross = @species * $number_of_chr;
    if ( $cross == @files ) {
        print "Perfect! All maf files ara pairwised.\n";
    }
    elsif ( $cross > @files ) {
        die "There are multiple species maf files. Not sure it works\n";
    }
    else {
        die "Please check maf files. It seemed there are redundances.\n";
    }
}

unless ($out_dir) {
    $out_dir = ucfirst $target_name . "vs" . uc roman( scalar @species );
}
unless ( -e $out_dir ) {
    mkdir $out_dir, 0777
        or die "Cannot create directory [$out_dir]: $!";
}

#----------------------------#
# matching mafs
#----------------------------#
#---
#Bur_0:
#  1: ~/data/alignment/arabidopsis19/AthvsBur_0/mafNet/chr1.net.maf
#  2: ...
#lyrata_65:
#  1: ~/data/alignment/arabidopsis19/AthvsLyrata_set01_4/mafNet/chr1.net.maf
#  2: ...
my $file_of = {};
my @chr_names;
{
    my @chr_string;
    for my $name (@species) {
        my @species_files;
        for my $key (%species_of) {
            my ($dummy) = grep { $_ ne $target_name } @{ $species_of{$key} };
            if ( $dummy and $name eq $dummy ) {
                push @species_files, $key;
            }
        }

        my @chrs = map { basename $_ , $suffix }
            @species_files;    # strip dir and suffix
        while (1) {
            my $lcss = lcss(@chrs);
            last unless $lcss;
            print "LCSS [$lcss]\n";
            my $rx = quotemeta $lcss;
            $chrs[$_] =~ s/$rx// for 0 .. $#chrs;
        }
        $file_of->{$name} = { zip( @chrs, @species_files ) };
        push @chr_string, join " ", sort @chrs;
    }

    if ( scalar uniq(@chr_string) == 1 ) {
        @chr_names = split / /, $chr_string[0];
        if ( @chr_names == $number_of_chr ) {
            print "Check maf filenames OK\n";
            print "@chr_names\n";
        }
        else {
            print "chr_string\t@chr_string\n";
            print "number_of_chr\t$number_of_chr\n";
            die "Check maf filenames FAILED\n";
        }
    }
    else {
        print "@chr_string\n";
        die "Check maf filenames FAILED\n";
    }
}

#----------------------------#
# Finally, multiz comes
#----------------------------#

my $worker = sub {
    my $job = shift;
    my $opt = shift;

    my $chr_name = $job;

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
            $maf1     = $file_of->{$species1}{$chr_name};
        }
        else {
            $maf1 = $maf_step;
        }
        if ( !defined $species2 ) {
            $species2 = shift @species_copy;
            $maf2     = $file_of->{$species2}{$chr_name};
        }

        my $out1 = "$out_dir/$chr_name.out1";
        my $out2 = "$out_dir/$chr_name.out2";

        $maf_step
            = @species_copy
            ? "$out_dir/chr$chr_name.step$step$suffix"
            : "$out_dir/chr$chr_name$suffix";

        # here we set out1 and out2 to discard unused synteny
        # Omit out1 and out2, unused synteny will be printed to stdout and
        # reused by following multiz processes
        print "Run multiz...\n";
        my $cmd
            = "$multiz_bin/multiz" 
            . " $maf1" 
            . " $maf2" . " 1 " 
            . " $out1"
            . " $out2"
            . " > $maf_step";
        exec_cmd($cmd);
        print "Step [$step] .maf file generated.\n\n";

        $str .= "chr$chr_name.step$step,";
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

    print "Clean temp files.\n";
    remove("$out_dir/$chr_name.out1");
    remove("$out_dir/$chr_name.out2");
    remove("$out_dir/chr$chr_name.step*");
    
    my $cmd = "gzip " . "$out_dir/chr$chr_name$suffix";
    exec_cmd($cmd);

    print $str;
    open my $out_fh, ">", "$out_dir/chr$chr_name.temp.csv";
    print {$out_fh} $str;
    close $out_fh;

    return;
};

my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@chr_names,
    code     => $worker,
);
$run->run;

{    # summary
    my $cmd = "echo -e 'step,spe1,spe2,maf1,maf2,out1,out2,size,per_size'"
        . " > $out_dir/steps.csv";
    exec_cmd($cmd);

    $cmd
        = "find $out_dir -type f -name '*.temp.csv'"
        . " | sort -n"
        . " | xargs cat >> $out_dir/steps.csv";
    exec_cmd($cmd);

    $cmd = "find $out_dir -type f -name '*.temp.csv'" . " | xargs rm";
    exec_cmd($cmd);
}

#----------------------------#
# compress
#----------------------------#
if ($gzip) {
    print "Restore compressed state of maf files.\n";
    print "Compress...\n";
    
    my $worker = sub {
        my $job = shift;

        my $file = $job;
        my $cmd  = "gzip $file";
        exec_cmd($cmd);

        return;
    };

    my @jobs = sort @files;
    my $run  = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;
}

print "\n";
print "Runtime ", duration( time - $start_time ), ".\n";
print "=" x 30, "\n";

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

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
            $distance_of->{ $node->get_name }
                = $node->calc_patristic_distance($target);
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

# comes from
# http://stackoverflow.com/questions/499967/how-do-i-determine-the-longest-similar-portion-of-several-strings
sub lcss {
    return '' unless @_;
    return $_[0] if @_ == 1;
    my $i          = 0;
    my $first      = shift;
    my $min_length = length($first);
    for (@_) {
        $min_length = length($_) if length($_) < $min_length;
    }
INDEX: for my $ch ( split //, $first ) {
        last INDEX unless $i < $min_length;
        for my $string (@_) {
            last INDEX if substr( $string, $i, 1 ) ne $ch;
        }
    }
    continue { $i++ }
    return substr $first, 0, $i;
}

__END__

=head1 NAME

    mz.pl - Multiz step by step

=head1 SYNOPSIS

    mz.pl -dt <one target dir or file> -dq <one query dir or file> [options]
      Options:
        -?, --help              brief help message
        --man                   full documentation

      Run in parallel mode
        -p, --paralle           number of child processes

      Fasta dirs  
        -dt, --dir_target       dir of target fasta files
        -dq, --dir_query        dir of query fasta files

      Output .lav and .axt
        -dl, --dir_lav          where .lav and .axt files storess
    
        perl mz.pl  -d ~/data/alignment/arabidopsis19/AthvsLyrata/ \
                    -d ~/data/alignment/arabidopsis19/AthvsBur_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsZu_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsNo_0/ \
                    -d ~/data/alignment/arabidopsis19/AthvsLer_0/ \
                    --tree ~/data/alignment/arabidopsis19/19way.nwk \
                    --parallel 8 \
                    -syn

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION


=cut
