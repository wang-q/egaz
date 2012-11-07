#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server     = $Config->{database}{server};
my $port       = $Config->{database}{port};
my $username   = $Config->{database}{username};
my $password   = $Config->{database}{password};
my $db         = $Config->{database}{db};
my $ensembl_db = $Config->{database}{ensembl};

# write_axt parameter
my $length_threshold = $Config->{write}{feature_threshold};
my $feature          = $Config->{write}{feature};

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    's|server=s'    => \$server,
    'P|port=i'      => \$port,
    'd|db=s'        => \$db,
    'u|username=s'  => \$username,
    'p|password=s'  => \$password,
    'e|ensembl=s'   => \$ensembl_db,
    'l|lt|length=i' => \$length_threshold,
    'feature=s'     => \$feature,
    'parallel=i'    => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write .axt files from $db...");

my @jobs;
{    # create alignDB object for this scope
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # select all target chromosomes in this database
    my @chrs = @{ $obj->get_chrs('target') };
    @jobs = @chrs;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job = shift;
    my $chr = $job;
    local $| = 1;

    #----------------------------#
    # Init objects
    #----------------------------#
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    # Database handler
    my $dbh = $obj->dbh;

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    # position finder
    my $pos_obj = AlignDB::Position->new( dbh => $dbh );

    # get target and query names
    my ( $target_name, $query_name ) = $obj->get_names;

    my ( $chr_id, $chr_name, $chr_length ) = @{$chr};
    print Dump {
        chr_id     => $chr_id,
        chr_name   => $chr_name,
        chr_length => $chr_length
    };

    # prepare axt summary line
    my $serial = 0;
    my $dir    = "${db}_${feature}";
    mkdir $dir, 0777 unless -e $dir;
    my $filename = "$dir/$chr_name.$feature.axt";

    # for each align
    my @align_ids = @{ $obj->get_align_ids_of_chr($chr_id) };
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};
        my $target_runlist   = $target_info->{seq_runlist};

        # query
        my $query_info      = $obj->get_query_info($align_id);
        my $query_chr_name  = $query_info->{chr_name};
        my $query_chr_start = $query_info->{chr_start};
        my $query_chr_end   = $query_info->{chr_end};
        my $query_runlist   = $query_info->{seq_runlist};
        my $query_strand    = $query_info->{query_strand};

        my ( $target_seq, $query_seq ) = @{ $obj->get_seqs($align_id) };

        # make a new ensembl slice object
        my $ensembl_chr_name = $target_chr_name;
        $ensembl_chr_name =~ s/chr0?//i;

        #print "ensembl_chr_name $ensembl_chr_name\n";
        eval {
            $ensembl->set_slice( $ensembl_chr_name, $target_chr_start,
                $target_chr_end );
        };
        if ($@) {
            warn "Can't get annotation\n";
            next;
        }

        my $slice       = $ensembl->slice;
        my $ftr_chr_set = $slice->{"_$feature\_set"};

        if ( $ftr_chr_set->empty ) {
            print "  No $feature, jump to next alignment\n";
            next;
        }

        my $target_set = AlignDB::IntSpan->new($target_runlist);
        my $query_set  = AlignDB::IntSpan->new($query_runlist);

        my $ftr_chr_runlist = $ftr_chr_set->run_list;
        my @ftr_chrs = split /(\,|-)/, $ftr_chr_runlist;
        my @ftr_aligns
            = map { /^\d+$/ ? $pos_obj->at_align( $align_id, $_ ) : $_; }
            @ftr_chrs;
        my $ftr_align_set = AlignDB::IntSpan->new( join '', @ftr_aligns );

        my @ftr_segments = split ",", $ftr_align_set->run_list;
        foreach (@ftr_segments) {
            next if /^\d+$/;
            my ( $seg_start, $seg_end ) = split "-";
            my $seg_length = $seg_end - $seg_start + 1;
            next if $seg_length < $length_threshold;

            # align coordinates to target & query chromosome coordinates
            my $target_seg_start
                = $pos_obj->at_target_chr( $align_id, $seg_start );
            my $target_seg_end = $pos_obj->at_target_chr( $align_id, $seg_end );
            my $query_seg_start
                = $pos_obj->at_query_chr( $align_id, $seg_start );
            my $query_seg_end = $pos_obj->at_query_chr( $align_id, $seg_end );
            if ( $query_strand eq '-' ) {
                ( $query_seg_start, $query_seg_end )
                    = ( $query_seg_end, $query_seg_start );
            }
            my $score = $seg_length * 100;    # fake score

            # append axt file
            {
                open my $outfh, '>>', $filename;
                print {$outfh} "$serial";
                print {$outfh} " $target_chr_name";
                print {$outfh} " $target_seg_start $target_seg_end";
                print {$outfh} " $query_chr_name";
                print {$outfh} " $query_seg_start $query_seg_end";
                print {$outfh} " $query_strand $score\n";
                print {$outfh}
                    substr( $target_seq, $seg_start - 1, $seg_length ), "\n";
                print {$outfh}
                    substr( $query_seq, $seg_start - 1, $seg_length ), "\n";
                print {$outfh} "\n";
                close $outfh;
            }

            $serial++;
        }

        print "  finish write axt file\n";
    }
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;

__END__

=head1 NAME

    write_axt.pl - extract sequence of a certain feature from alignDB

=head1 SYNOPSIS

    write_axt.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --server            MySQL server IP/Domain name
        --db                database name
        --username          username
        --password          password
        --ensembl           ensembl database name

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

