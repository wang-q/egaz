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
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

my $query_dir        = '';
my $length_threshold = 5000;
my $wrap_length      = 50;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=s'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'query_dir=s'  => \$query_dir,
    'length=i'     => \$length_threshold,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

pod2usage( -exitstatus => 0, -verbose => 2 ) unless $query_dir;

my @target_chr_ids = (59);
my @query_chr_ids  = (209);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Write unmatched query sequence from $db...");

# alignDB object
my $alignDB_obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

my $dbh = $alignDB_obj->dbh();

# query's chromosomal location
#
my $query_query = q{
    SELECT  s.chr_start,
            s.chr_end
    FROM    query q, sequence s,
            (SELECT a.align_id
            FROM align a, target t, sequence s
            WHERE a.align_id = t.align_id
            AND t.seq_id = s.seq_id
            AND s.chr_id = ?
            ORDER BY a.align_id) a
    WHERE q.seq_id = s.seq_id
    AND s.chr_id = ?
    AND q.align_id = a.align_id
};
my $query_sth = $dbh->prepare($query_query);

# query chromosomal info
#
my $query_chr_query = q{
    SELECT  c.chr_name,
            c.chr_length
    FROM    chromosome c
    WHERE c.chr_id = ?
};
my $query_chr_sth = $dbh->prepare($query_chr_query);

# store all query chr locations
my $query_set = {};

foreach my $target_chr_id (@target_chr_ids) {
    foreach my $query_chr_id (@query_chr_ids) {

        $query_chr_sth->execute($query_chr_id);
        my ( $chr_name, $chr_length ) = $query_chr_sth->fetchrow_array;
        print "Processing query $chr_name\n";

        my $chr_seq;

        open CHRFH, "<$query_dir/$chr_name.fa";
        my $dummy = <CHRFH>;
        {
            local $/;
            $chr_seq = <CHRFH>;
        }
        $chr_seq =~ tr{/\\\r\n\b\f. }{}d;

        $chr_seq = uc $chr_seq;

        #open TEMP, ">temp.fasta";
        #print TEMP $chr_seq;
        #close TEMP;

        print "Read in chromosome $chr_name sequences\n";

        my $full_chr_set       = AlignDB::IntSpan->new("1-$chr_length");
        my $query_chr_set      = AlignDB::IntSpan->new();
        my $complement_chr_set = AlignDB::IntSpan->new();

        $query_sth->execute( $target_chr_id, $query_chr_id );
        while ( my @row = $query_sth->fetchrow_array ) {
            my ( $chr_start, $chr_end ) = @row;
            $query_chr_set->add("$chr_start-$chr_end");
        }

        $complement_chr_set = $full_chr_set->diff($query_chr_set);

        print "complement_chr_set generated\n";

        open OUTFH, ">$chr_name.complement.fa";
        my @spans = $complement_chr_set->spans();
        my $i     = 0;
    SEG: foreach my $span (@spans) {
            my $seg_start  = $span->[0];
            my $seg_end    = $span->[1];
            my $seg_length = $seg_end - $seg_start + 1;
            next if ( $seg_length < $length_threshold );

            my $seq = substr( $chr_seq, $seg_start - 1, $seg_length );
            my @sub_seqs = split /[Nn]{50,}/, $seq;

            foreach my $sub_seq (@sub_seqs) {
                my $sub_seq_length = length $sub_seq;
                next if ( $sub_seq_length < $length_threshold );
                my $n_number = $sub_seq =~ tr/Nn/Nn/;
                if ( $n_number / $sub_seq_length > 0.1 ) {
                    next;
                }
                my $header = ">$chr_name" . "_" . "$i";
                print OUTFH "$header\n";
                if ($wrap_length) {
                    for (
                        my $pos = 0;
                        $pos < $sub_seq_length;
                        $pos += $wrap_length
                        )
                    {
                        print OUTFH substr( $sub_seq, $pos, $wrap_length ),
                            "\n";
                    }
                }
                else {
                    print OUTFH "$seq\n";
                }

                $i++;
            }
        }
        close OUTFH;
    }
}

$stopwatch->end_message();
exit;

__END__


=head1 NAME

    alignDB_graph.pl - Generate graph for chromosome coverage in alignDB

=head1 SYNOPSIS

    alignDB_graph.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --class           GD or GD::SVG
       

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
