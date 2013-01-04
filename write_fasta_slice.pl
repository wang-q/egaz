#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Spec;
use File::Basename;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# write_axt parameter
my $yaml_dir = '.';

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'y|yaml_dir=s' => \$yaml_dir,
    'parallel=i'   => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write .fas files from $db...");

#----------------------------------------------------------#
# Write .axt files from alignDB
#----------------------------------------------------------#
my @yaml_files
    = File::Find::Rule->file->name( '*.yaml', '*.yml' )->in($yaml_dir);
printf "\n----Total YAML Files: %4s----\n\n", scalar @yaml_files;

for my $yaml_file ( sort @yaml_files ) {
    print "Loading $yaml_file\n";
    my ( $base, $dir ) = fileparse( $yaml_file, ".yaml", ".yml" );

    if ( $base =~ /^(chr\w+)/ ) {
        my $chr_name  = $1;
        my $runlist   = LoadFile($yaml_file);
        my $slice_set = AlignDB::IntSpan->new($runlist);
        my $outfile   = File::Spec->catfile( $dir, "$base.fas" );
        print "Write $outfile\n";
        write_slice( $chr_name, $slice_set, $outfile );
    }
    else {
        my $slice_set_of = LoadFile($yaml_file);

        my $worker = sub {
            my $chr_name  = shift;
            my $slice_set = AlignDB::IntSpan->new( $slice_set_of->{$chr_name} );
            my $outfile   = File::Spec->catfile( $dir, "$chr_name.fas" );
            write_slice( $chr_name, $slice_set, $outfile );
        };

        my $run = AlignDB::Run->new(
            parallel => $parallel,
            jobs     => [ sort keys %{$slice_set_of} ],
            code     => $worker,
        );
        $run->run;
    }
}

sub write_slice {
    my $chr_name  = shift;
    my $slice_set = shift;
    my $outfile   = shift;

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # position finder
    my $pos_obj = AlignDB::Position->new( dbh => $dbh );

    print "Write slice from $chr_name\n";
    print "Output file is $outfile\n";

    # alignment
    my @align_ids = @{ $obj->get_align_ids_of_chr_name($chr_name) };

    for my $align_id (@align_ids) {
        local $| = 1;
        print "Processing align_id $align_id\n";

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};
        my $target_runlist   = $target_info->{seq_runlist};

        my ( $target_seq, @query_seqs ) = @{ $obj->get_seqs($align_id) };
        my $ref_seq = $obj->get_seq_ref($align_id);

        my $target_set = AlignDB::IntSpan->new($target_runlist);

        my $align_chr_set
            = AlignDB::IntSpan->new("$target_chr_start-$target_chr_end");
        my $iset = $slice_set->intersect($align_chr_set);
        next if $iset->is_empty;

        # there may be two or more subslice intersect this alignment
        for my $ss_set ( $iset->sets ) {

            # rhs position set
            my $ss_start = $pos_obj->at_align( $align_id, $ss_set->min );
            my $ss_end   = $pos_obj->at_align( $align_id, $ss_set->max );
            next if $ss_start >= $ss_end;
            $ss_set = AlignDB::IntSpan->new("$ss_start-$ss_end");
            $ss_set = $ss_set->intersect($target_set);

            next if $ss_set->count <= 1;
            my ( $seg_start, $seg_end ) = ( $ss_set->min, $ss_set->max );
            my $seg_length = $seg_end - $seg_start + 1;

            # align coordinates to target & query chromosome coordinates
            my $target_seg_start
                = $pos_obj->at_target_chr( $align_id, $seg_start );
            my $target_seg_end = $pos_obj->at_target_chr( $align_id, $seg_end );

            # append fas file
            # S288C.chrI(+):27520-29557|species=S288C
            {
                print "Write fasta files: "
                    . "$target_chr_name:$target_seg_start-$target_seg_end"
                    . "\n";
                open my $outfh, '>>', $outfile;
                print {$outfh} ">ref\n";
                print {$outfh} substr( $ref_seq, $seg_start - 1, $seg_length ),
                    "\n";
                print {$outfh}
                    ">target.$target_chr_name(+):$target_chr_start-$target_chr_end|species=target\n";
                print {$outfh}
                    substr( $target_seq, $seg_start - 1, $seg_length ), "\n";

                for my $i ( 0 .. $#query_seqs ) {
                    print {$outfh} ">query$i\n";
                    print {$outfh}
                        substr( $query_seqs[$i], $seg_start - 1, $seg_length ),
                        "\n";
                }
                print {$outfh} "\n";
                close $outfh;
            }
        }
        print "  finish write fasta file\n";
    }
}

$stopwatch->end_message;

__END__

=head1 NAME

    write_fasta_slice.pl - extract alignment slices from MalignDB

=head1 SYNOPSIS

    write_fasta_slice.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        -y, --yaml_dir      dir of yaml

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

format 1: chr1.yml
  --- 1-25744,815056-817137
  
  output axt: chr1/chr1.axt
format 2: bin1.yml
  ---
  chr1: '1-25744,815056-817137'

  output axt: bin1/chr1.axt


=cut

