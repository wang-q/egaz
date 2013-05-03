#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use File::Spec;

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ensembl;

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

# removing $n integers from each end of each span of $set.
# If $n is negative, then -$n integers are added to each end of each span.
# use AlignDB::IntSpan->inset()
my $inset;

# write inverted sets
my $invert;

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'        => \$help,
    'man'           => \$man,
    's|server=s'    => \$server,
    'P|port=i'      => \$port,
    'u|username=s'  => \$username,
    'p|password=s'  => \$password,
    'd|db=s'        => \$db,
    'e|ensembl=s'   => \$ensembl_db,
    'l|lt|length=i' => \$length_threshold,
    'feature=s'     => \$feature,
    'inset=i'       => \$inset,
    'invert'        => \$invert,
    'parallel=i'    => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write slice files from $db...");

# output dir
my $dir = "${db}_${feature}";
$dir .= "_inset$inset" if $inset;
$dir .= "_invert"      if $invert;
mkdir $dir, 0777 unless -e $dir;

#----------------------------#
# Find all align_ids
#----------------------------#
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

    # ensembl handler
    my $ensembl = AlignDB::Ensembl->new(
        server => $server,
        db     => $ensembl_db,
        user   => $username,
        passwd => $password,
    );

    my ( $chr_id, $chr_name, $chr_length ) = @{$chr};
    print "id => $chr_id, name => $chr_name, length => $chr_length\n";

    # for each align
    my @align_ids     = @{ $obj->get_align_ids_of_chr($chr_id) };
    my $chr_ftr_set   = AlignDB::IntSpan->new;
    my $chr_align_set = AlignDB::IntSpan->new;
    for my $align_id (@align_ids) {
        $obj->process_message($align_id);

        # target
        my $target_info      = $obj->get_target_info($align_id);
        my $target_chr_name  = $target_info->{chr_name};
        my $target_chr_start = $target_info->{chr_start};
        my $target_chr_end   = $target_info->{chr_end};

        $chr_align_set->add("$target_chr_start-$target_chr_end");

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

        next unless $ftr_chr_set;
        next if $ftr_chr_set->is_empty;
        if ($inset) {
            $ftr_chr_set->inset($inset);
        }

        for my $set ( $ftr_chr_set->sets ) {
            next if $set->size < $length_threshold;
            $chr_ftr_set->add($set);
        }
    }

    # $chr_ftr_set should be subset of $chr_align_set
    $chr_ftr_set = $chr_ftr_set->intersect($chr_align_set);

    # the inverted set
    if ($invert) {
        $chr_ftr_set = $chr_align_set->diff($chr_ftr_set);
    }

    my $filename = File::Spec->catfile( $dir, "$chr_name.$feature.yml" );
    if ( $chr_ftr_set->is_not_empty ) {
        DumpFile( $filename, $chr_ftr_set->runlist );
        print "Finish merge $feature of $chr_name\n";
    }
    else {
        print "Finish $chr_name has no $feature\n";
    }
};

#----------------------------------------------------------#
# start
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

    write_runlist_feature.pl - extract runlists of a certain feature from alignDB

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

=cut

