#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Paralog;
use AlignDB::Stopwatch;

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

my $datalib = $Config->{paralog}{datalib};

my $man  = 0;
my $help = 0;

$|++;

GetOptions(
    'help|?'       => \$help,
    'man|m'        => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'da|datalib=s' => \$datalib,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Write paralog alignment files of $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# build $paralog_obj with parameters from alignDB.ini
my $paralog_obj = AlignDB::Paralog->new;

# we need sequence, so couldn't use megablast
$paralog_obj->use_megablast(0);

# Normal blast result, pairwise alignment
$paralog_obj->alignment_view(0);

# the first hit is the sequence itself, so we want two hits
$paralog_obj->alignment_limit(2);
$paralog_obj->description_limit(2);

#----------------------------------------------------------#
# write_full_paralog alignment files
#----------------------------------------------------------#
# alignments
my $align_query = q{
    SELECT a.align_id, a.align_paralog
    FROM align a
    WHERE a.align_paralog * a.align_length > 10000
};
my $align_sth = $dbh->prepare($align_query);

# sequence
my $seq_query = q{
    SELECT  c.chr_name,
            s.chr_start,
            s.chr_end,
            s.seq_length,
            t.target_seq,
            t.target_runlist
    FROM target t, sequence s, chromosome c
    WHERE t.seq_id = s.seq_id
    AND s.chr_id = c.chr_id
    AND t.align_id = ?
};
my $seq_sth = $dbh->prepare($seq_query);

$align_sth->execute;
my $axt_serial = 0;

# for each align
while ( my @row = $align_sth->fetchrow_array ) {
    my ( $align_id, $align_paralog ) = @row;
    print "\nProcessing align_id $align_id,", " align_paralog $align_paralog\n";

    $seq_sth->execute($align_id);

    my ($target_chr_name, $target_chr_start, $target_chr_end,
        $target_length,   $target_seq,       $target_runlist
    ) = $seq_sth->fetchrow_array;
    $target_seq =~ tr/-//d;

    my ( $longest_align, $longest_hit_name, $longest_hit_strand );
    eval {

        # "out of memory" error
        my $report_ref = $paralog_obj->web_blast( \$target_seq );
        ( $longest_align, $longest_hit_name, $longest_hit_strand )
            = $paralog_obj->longest_paralog( $report_ref, $Config );

        # "undefined scalar as ref" error
        # Does a s/$arg1/$arg2/ on the sequences. Useful for gap characters
        $longest_align->map_chars( '\.', '-' );
        $longest_align->uppercase;
    };
    if ($@) {
        print " " x 2, "Blast process error, jump to next\n";
        next;
    }

    my ($query_obj) = $longest_align->each_seq_with_id("query");
    my ($hit_obj)   = $longest_align->each_seq_with_id($longest_hit_name);

    #print Dump([$query_header, $query_obj,$hit_header,$hit_obj ]);

    # output a fasta alignment file
    {
        print " " x 2, "Write FASTA file\n";

        my $query_header
            = "align_id_$align_id|" . $query_obj->start . "-" . $query_obj->end;
        my $hit_header
            = "$longest_hit_name|" . $hit_obj->start . "-" . $hit_obj->end;

        unless ( -e $db ) {
            mkdir $db, 0777
                or die "Cannot create \"$db\" directory: $!";
        }
        my $outfile = "$db/" . $align_id . ".fasta";
        open my $fasta_fh, '>', $outfile
            or die("Cannot open OUT file $outfile");
        print {$fasta_fh} ">$query_header\n";
        print {$fasta_fh} $query_obj->seq, "\n";
        print {$fasta_fh} ">$hit_header\n";
        print {$fasta_fh} $hit_obj->seq, "\n";
        close $fasta_fh;

        print " " x 4, "$outfile\n";
    }

    # output an axt alignment file
    {
        print " " x 2, "Write AXT file\n";

        my $query_chr_name = $target_chr_name;
        my $query_start    = $target_chr_start + $query_obj->start - 1;
        my $query_end      = $target_chr_start + $query_obj->end - 1;
        my $score = ( $query_end - $query_start + 1 ) * 100;    # sham score

        my $hit_chr_name = $longest_hit_name;
        my $hit_start    = $hit_obj->start;
        my $hit_end      = $hit_obj->end;

        #$hit_obj doesn't contain strand info
        my $hit_strand = $longest_hit_strand;

        #print Dump {
        #    hit_obj  => $hit_obj,
        #    hit_name => $hit_obj->id,
        #};

        unless ( -e $db ) {
            mkdir $db, 0777
                or die "Cannot create \"$db\" directory: $!";
        }
        my $outfile = "$db/" . $db . ".paralog.axt";

        # append axt file
        open my $axt_fh, '>>', $outfile
            or die "Cannot open OUT file $outfile";

        print {$axt_fh} "$axt_serial";
        print {$axt_fh} " $query_chr_name $query_start $query_end";
        print {$axt_fh} " $hit_chr_name $hit_start $hit_end";
        print {$axt_fh} " $hit_strand $score\n";
        print {$axt_fh} $query_obj->seq, "\n";
        print {$axt_fh} $hit_obj->seq, "\n";
        print {$axt_fh} "\n";
        close $axt_fh;

        print " " x 4, "$outfile\n";

        $axt_serial++;
    }
}

$align_sth->finish;

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    write_full_paralog.pl - After update_align_paralog.pl, write
                            longest_paralog-target alignments to files
    
=head1 SYNOPSIS

    write_full_paralog.pl [options]
     Options:
       --help               brief help message
       --man                full documentation
       --server             MySQL server IP/Domain name
       --db                 database name
       --username           username
       --password           password
       --datalib|da         blast database
       
    write_full_paralog.pl -d=Nipvs9311 -da=nip_chro

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

