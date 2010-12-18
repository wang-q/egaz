#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Run;
use AlignDB::Util qw(:all);

use Time::Duration;
use File::Find::Rule;
use File::Basename;
use Path::Class;
use String::Compare;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# executable file location
my $path_blastz    = "$FindBin::Bin/blastz";
my $path_normalize = "perl $FindBin::Bin/normalize_lav.pl";
my $path_lav2axt   = "perl $FindBin::Bin/lav2axt.pl";

# run in parallel mode
my $parallel = 1;

# use a specified set to run blastz
my $specified;

# use parted seqs, default is false
my ( $t_parted, $q_parted );

# relationship of target and query is one to one
my $paired;

# dirs
my ( $dir_target, $dir_query );
my $dir_lav = ".";

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

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'              => \$help,
    'man'                 => \$man,
    'pb|path_blastz=s'    => \$path_blastz,
    'pn|path_normalize=s' => \$path_normalize,
    'pl|path_lav2axt=s'   => \$path_lav2axt,
    'dt|dir_target=s'     => \$dir_target,
    'dq|dir_query=s'      => \$dir_query,
    'dl|dir_lav=s'        => \$dir_lav,
    's|specified=s'       => \$specified,
    'p|parallel=i'        => \$parallel,
    'tp|t_parted'         => \$t_parted,
    'qp|q_parted'         => \$q_parted,
    'paired'              => \$paired,
    'O=s'                 => \$opt{O},
    'E=s'                 => \$opt{E},
    'Q=s'                 => \$opt{Q},
    'C=s'                 => \$opt{C},
    'T=s'                 => \$opt{T},
    'M=s'                 => \$opt{M},
    'K=s'                 => \$opt{K},
    'L=s'                 => \$opt{L},
    'H=s'                 => \$opt{H},
    'Y=s'                 => \$opt{Y},
    'Z=s'                 => \$opt{Z},
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
# predefined parameter sets
# if given a specified set, run with this parameter set
my $parameters = LoadFile("$FindBin::Bin/parameters.yaml");
if ($specified) {
    print "*** Use parameter set $specified\n";
    my $para_set = $parameters->{$specified};
    foreach my $key ( keys %{$para_set} ) {
        next if $key eq "comment";
        $opt{$key} = $para_set->{$key};
        print "$key = $opt{$key} ";
    }
    print "\n";
}

# matrix file
$opt{Q} = "$FindBin::Bin/matrix/" . $opt{Q} if $opt{Q};

# make lav dir
unless ( -e $dir_lav ) {
    mkdir $dir_lav, 0777
        or die "Cannot create \"$dir_lav\" directory: $!";
}

#----------------------------------------------------------#
# findfile section
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
# blastz section
#----------------------------------------------------------#
{
    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my ( $target, $query ) = split /\|/, $job;

        print "Run blastz...\n";

        # naming the .lav file
        my $t_base = basename($target);
        $t_base =~ s/\..+?$//;
        my $q_base = basename($query);
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

        my $bz_cmd = "$path_blastz $target $query";
        for my $key ( keys %{$opt} ) {
            my $value = $opt->{$key};
            if ( defined $value ) {
                $bz_cmd .= " $key=$value";
            }
        }
        $bz_cmd .= " > $lav_file";
        exec_cmd($bz_cmd);
        print ".lav file generated.\n\n";

        return;
    };

    # All jobs to be done
    my @jobs;
    if ($paired) {

        # use the most similar chr name
        for my $target_file ( sort @target_files ) {
            my $t_base = basename($target_file);
            my ($query_file)
                = map { $_->[0] }
                sort  { $b->[1] <=> $a->[1] }
                map { [ $_, compare( basename($_), $t_base ) ] } @query_files;
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

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
        opt      => \%opt,
    );
    $run->run;

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
        %t_length = read_length_csv( $dir_target, 'length.csv' );
    }
    if ($q_parted) {
        %q_length = read_length_csv( $dir_query, 'length.csv' );
    }

    for my $file (@lav_files) {
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
            = "$path_normalize"
            . " -0 $t_len -1 $q_len"
            . " -i $file -o $outfile";

        exec_cmd($cmd);
        unlink $file;
        print ".lav file normalized.\n\n";
    }
}

#----------------------------------------------------------#
# lav2axt section
#----------------------------------------------------------#
{
    my @lav_files = File::Find::Rule->file->name('*.lav')->in($dir_lav);
    printf "\n----%4s .lav files to be converted ----\n", scalar @lav_files;

    my $worker = sub {
        my $job = shift;
        my $opt = shift;

        my $file = $job;

        print "Run lav2axt...\n";
        my $cmd = "$path_lav2axt -l $file ";
        exec_cmd($cmd);
        print ".axt file generated.\n\n";

        return;
    };

    my @jobs = sort @lav_files;

    my $start_time = time;
    print "\n", "=" x 30, "\n";
    print "Processing...\n";

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@jobs,
        code     => $worker,
    );
    $run->run;

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

sub read_length_csv {
    my $fh = file(@_)->openr;
    my %length_of;
    while (<$fh>) {
        chomp;
        my ( $key, $value ) = split /,/;
        $length_of{$key} = $value;
    }

    return %length_of;
}

__END__

=head1 NAME

    bz.pl - execute blastz and lav2axt against two directories
    Fasta file naming rules: "seqfile[from,to]"
    Lav file naming rules: "[target]vs[query].N.lav"

=head1 SYNOPSIS

    bz.pl -dt <one target dir or file> -dq <one query dir or file> [options]
      Options:
        -?, --help              brief help message
        --man                   full documentation

      Run in parallel mode
        -p, --paralle           number of child processes

      Fasta dirs  
        -dt, --dir_target       dir of target fasta files
        -dq, --dir_query        dir of query fasta files

      Output .lav and .axt
        -dl, --dir_lav          where .lav and .axt files stores

      Executable files:
        -pb, --path_blastz      path to blastz executable file
        -pn, --path_normalize   path to normalize_lav.pl
        -pl, --path_lav2axt     path to lav2axt.pl

      Predefined parameter set:  
        -s,  --specified        use a predefined parameter set
        
      Relationship
        --paired                relationship of target and query is one to one

      Scoring parameters:
        -O                      gap-open penalty
        -E                      gap-extension penalty
        -Q                      matrix file

      Aligning parameters: 
        -C                      chain option
        -T                      words option
        -M                      mask any base in seq1 hit this many times

      Droping hsp parameters:
        -K                      threshold for MSPs for the first pass
        -L                      threshold for gapped alignments
                                for the second pass
        -H                      threshold to be interpolated between
                                alignments
        -Y                      X-drop parameter for gapped extension

      Speedup parameters:
        -Z                      increment between successive words
    
    >perl part_seq.pl -in t\RM11 -out t\parted -chunk 500000
    >dos2unix *.fa
    >perl bz.pl -dt t\S288C\chr01.fa -dq t\parted -s set01 -dl test -qp
    
    >perl part_seq.pl -in t\S288C -out t\S288C_parted -chunk 500000
    >perl bz.pl -dt t\S288C_parted -dq t\RM11\rm11.fa -s set01 -dl test -tp -p 4

=head1 OPTIONS

=over 4

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

Blastz will take the first sequence in target fasta file and all sequences in
query fasta file.

So if there are mutiple query files, the program uses the largest one. And all
target files in dir_target will be processed.

So, if there are combined fasta files and multi fasta files coexisting in the
target directory, just delete the axt file matched with the combined fasta
filename.

=cut
