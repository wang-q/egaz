#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use App::RL::Common;
use File::Find::Rule;
use MCE;
use Path::Tiny;
use String::Compare;
use Time::Duration;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

bz.pl - execute lastz and lav2axt against two paths

=head1 SYNOPSIS

    perl bz.pl -dt <one target dir or file> -dq <one query dir or file> [options]
      Options:
        --help          -?          brief help message

      Running mode
        --parallel      -p  INT     run in parallel mode, [1]
        --noaxt                     don't generate .axt files

      Fasta dirs  
        --dir_target    -dt STR     dir of target fasta files
        --dir_query     -dq STR     dir of query fasta files

      Output .lav and .axt
        --dir_lav       -dl STR     where .lav and .axt files stores

      Inputs are parted or not
        --t_parted      -tp         use parted seqs, default is false
        --q_parted  `   -qp

      Predefined parameter set:  
        --specified     -s  STR     use a predefined parameter set
        
      Relationship
        --paired                    relationship of target and query is one to one
        --is_self                   self-alignment
                                    http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.01.50/README.lastz-1.01.50.html#ex_self

      Scoring parameters:
        -O                  INT     gap-open penalty
        -E                  INT     gap-extension penalty
        -Q                  INT     matrix file

      Aligning parameters: 
        -C                  INT     chain option
        -T                  INT     words option
        -M                  INT     mask any base in seq1 hit this many times

      Droping hsp parameters:
        -K                  INT     threshold for MSPs for the first pass
        -L                  INT     threshold for gapped alignments for the second pass
        -H                  INT     threshold to be interpolated between alignments
        -Y                  INT     X-drop parameter for gapped extension

      Speedup parameters:
        -Z                  INT     increment between successive words

    perl part_seq.pl -i t/S288C -o t/S288C_parted -chunk 500000
    perl bz.pl -dt t/S288C_parted -dq t/RM11/RM11.fa -s set01 -dl t/S288CvsRM11_df_tp -tp -p 1

=head1 DESCRIPTION

Lastz will take the first sequence in target fasta file and all sequences in query fasta file.

So if there are mutiple query files, the program uses the largest one. And all target files in
dir_target will be processed.

So, if there are combined fasta files and multi fasta files coexisting in the target directory, just
delete the axt file matched with the combined fasta filename.

Fasta file naming rules: "seqfile[from,to]"

Lav file naming rules: "[target]vs[query].N.lav"

=cut

my %opt = (
    O => undef,
    E => undef,
    Q => undef,
    C => undef,
    T => undef,
    M => undef,
    K => undef,
    L => undef,
    H => undef,
    Y => undef,
    Z => undef,
);

GetOptions(
    'help|?'          => sub { Getopt::Long::HelpMessage(0) },
    'dir_target|dt=s' => \my $dir_target,
    'dir_query|dq=s'  => \my $dir_query,
    'dir_lav|dl=s' => \( my $dir_lav  = '.' ),
    'parallel|p=i' => \( my $parallel = 1 ),
    'specified|s=s' => \my $specified,
    't_parted|tp'   => \my $t_parted,
    'q_parted|qp'   => \my $q_parted,
    'paired'        => \my $paired,
    'is_self'       => \my $is_self,
    'noaxt'         => \my $noaxt,
    'O=s'           => \$opt{O},
    'E=s'           => \$opt{E},
    'Q=s'           => \$opt{Q},
    'C=s'           => \$opt{C},
    'T=s'           => \$opt{T},
    'M=s'           => \$opt{M},
    'K=s'           => \$opt{K},
    'L=s'           => \$opt{L},
    'H=s'           => \$opt{H},
    'Y=s'           => \$opt{Y},
    'Z=s'           => \$opt{Z},
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
# predefined parameter sets
# if given a specified set, run with this parameter set
my $parameters = LoadFile("$FindBin::RealBin/parameters.yaml");
if ($specified) {
    print "*** Use parameter set $specified\n";
    my $para_set = $parameters->{$specified};
    for my $key ( keys %{$para_set} ) {
        next if $key eq "comment";
        next if defined $opt{$key};
        $opt{$key} = $para_set->{$key};
        print "$key = $opt{$key} ";
    }
    print "\n";
}

# matrix file
$opt{Q} = "$FindBin::RealBin/matrix/" . $opt{Q} if $opt{Q};

# make lav dir
path($dir_lav)->mkpath;

#----------------------------------------------------------#
# find file section
#----------------------------------------------------------#
my ( @target_files, @query_files );
if ($dir_target) {
    if ($t_parted) {
        @target_files = File::Find::Rule->file->name('*[*]')->in($dir_target);
    }
    else {
        @target_files = File::Find::Rule->file->name('*.fa')->in($dir_target);
    }
}
else {
    die "You should provide a dir or file for target.\n";
}

if ($dir_query) {
    if ($q_parted) {
        @query_files = File::Find::Rule->file->name('*[*]')->in($dir_query);
    }
    else {
        @query_files = File::Find::Rule->file->name('*.fa')->in($dir_query);
    }
}
else {
    die "You should provide a dir or file for query.\n";
}

printf "\n----%4s .fa files for target----\n", scalar @target_files;
printf "\n----%4s .fa files for query----\n",  scalar @query_files;

#----------------------------------------------------------#
# lastz section
#----------------------------------------------------------#
{
    my $worker = sub {
        my ( $self, $chunk_ref, $chunk_id ) = @_;

        my $job = $chunk_ref->[0];

        my ( $target, $query ) = split /\|/, $job;

        print "Run lastz...\n";

        # naming the .lav file
        # remove .fa or .fa[1,10000]
        my $t_base = path($target)->basename;
        $t_base =~ s/\..+?$//;
        my $q_base = path($query)->basename;
        $q_base =~ s/\..+?$//;

        my $lav_file;
        my $i = 0;
        while (1) {
            my $file = "[${t_base}]vs[${q_base}].$i.lav";
            $file = $dir_lav . "/" . $file;
            if ( !-e $file ) {
                $lav_file = $file;
                last;
            }
            $i++;
        }

        my $bz_cmd = "lastz $target $query";
        if ( $is_self and $target eq $query ) {
            $bz_cmd = "lastz $target --self";
        }

        for my $key ( keys %opt ) {
            my $value = $opt{$key};
            if ( defined $value ) {
                $bz_cmd .= " $key=$value";
            }
        }
        $bz_cmd .= " --ambiguous=iupac";
        $bz_cmd .= " > $lav_file";
        exec_cmd($bz_cmd);
        printf "\n.lav file generated. [%s]\n\n", path($lav_file)->basename;

        return;
    };

    # All jobs to be done
    my @jobs;
    if ($paired) {    # use the most similar chr name
        for my $target_file ( sort @target_files ) {
            my $t_base = path($target_file)->basename;
            my ($query_file) = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [ $_, compare( path($_)->basename, $t_base ) ] }
                @query_files;
            push @jobs, "$target_file|$query_file";
        }
    }
    else {
        for my $target_file ( sort @target_files ) {
            for my $query_file ( sort @query_files ) {
                push @jobs, "$target_file|$query_file";
            }
        }
    }

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach( \@jobs, $worker );

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

#----------------------------------------------------------#
# normalize section
#----------------------------------------------------------#
if ( $t_parted or $q_parted ) {
    my @lav_files = File::Find::Rule->file->name('*.lav')->in($dir_lav);
    printf "\n----%4s .lav files to be normalized ----\n", scalar @lav_files;

    my ( %t_length, %q_length );
    if ($t_parted) {
        %t_length = %{
            App::RL::Common::read_sizes(
                path( $dir_target, 'chr.sizes' )->stringify
            )
        };
    }
    if ($q_parted) {
        %q_length = %{
            App::RL::Common::read_sizes(
                path( $dir_query, 'chr.sizes' )->stringify
            )
        };
    }

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @lav_files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file = $chunk_ref->[0];

            $file =~ /\[(.+?)\]vs\[(.+?)\]/;
            my $t_name = $1;
            my $q_name = $2;

            if ( !defined $t_name or !defined $q_name ) {
                die ".lav naming error.\n";
            }

            my $outfile = $file;
            $outfile =~ s/\.lav$/\.norm\.lav/;

            my $t_len = $t_parted ? $t_length{$t_name} : 0;
            my $q_len = $q_parted ? $q_length{$q_name} : 0;

            print "Run normalize lav...\n";
            my $cmd
                = "perl $FindBin::RealBin/normalize_lav.pl"
                . " -0 $t_len -1 $q_len"
                . " -i $file -o $outfile";

            exec_cmd($cmd);
            path($file)->remove;
            print ".lav file normalized.\n\n";
        }
    );
}

#----------------------------------------------------------#
# lav2axt section
#----------------------------------------------------------#
if ( !$noaxt ) {
    my @lav_files = File::Find::Rule->file->name('*.lav')->in($dir_lav);
    printf "\n----%4s .lav files to be converted ----\n", scalar @lav_files;

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $mce = MCE->new( chunk_size => 1, max_workers => $parallel, );
    $mce->foreach(
        [ sort @lav_files ],
        sub {
            my ( $self, $chunk_ref, $chunk_id ) = @_;

            my $file = $chunk_ref->[0];

            print "Run lav2axt...\n";
            my $cmd = "perl $FindBin::RealBin/lav2axt.pl -i $file ";
            exec_cmd($cmd);
            print ".axt file generated.\n\n";
        }
    );

    print "\n";
    print "Runtime ", duration( time - $start_time ), ".\n";
    print "=" x 30, "\n";
}

exit;

sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

__END__
