package MyUtil;
use strict;
use warnings;
use autodie;

use Carp;
use Path::Tiny;
use Tie::IxHash;
use List::MoreUtils qw(minmax);

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            string_to_set set_to_string change_strand read_sizes revcom exec_cmd get_seq_faidx
            decode_header change_name_chopped
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub string_to_set {
    my $node = shift;

    my ( $chr, $runlist ) = split /:/, $node;
    my $strand = "+";
    if ( $chr =~ /\((.+)\)/ ) {
        $strand = $1;
        $chr =~ s/\(.+\)//;
    }
    my $set = AlignDB::IntSpan->new($runlist);

    return ( $chr, $set, $strand );
}

sub set_to_string {
    my $node = shift;

    my $string = $node->[0] . "(" . $node->[2] . "):" . $node->[1]->runlist;

    return $string;
}

sub change_strand {
    my $strand = shift;

    if ( $strand eq '+' ) {
        return '-';
    }
    elsif ( $strand eq '-' ) {
        return '+';
    }
    else {
        return $strand;
    }
}

sub read_sizes {
    my $file       = shift;
    my $remove_chr = shift;

    my @lines = path($file)->lines( { chomp => 1 } );
    tie my %length_of, "Tie::IxHash";
    for (@lines) {
        my ( $key, $value ) = split /\t/;
        $key =~ s/chr0?// if $remove_chr;
        $length_of{$key} = $value;
    }

    return \%length_of;
}

sub revcom {
    my $seq = shift;

    _ref2str( \$seq );
    $seq =~ tr/ACGTMRWSYKVHDBNacgtmrwsykvhdbn-/TGCAKYWSRMBDHVNtgcakywsrmbdhvn-/;
    my $seq_rc = reverse $seq;

    return $seq_rc;
}

# in situ convert reference of string to string
# For the sake of efficiency, the return value should be discarded
sub _ref2str {
    my $ref = shift;

    if ( ref $ref eq "REF" ) {
        $$ref = $$$ref;    # this is very weird, but it works
    }

    unless ( ref $ref eq "SCALAR" ) {
        carp "Wrong parameter passed\n";
    }

    return $ref;
}

sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd;
}

sub get_seq_faidx {
    my $genome   = shift;
    my $location = shift;    # I:1-100

    my $cmd = sprintf "samtools faidx %s %s", $genome, $location;
    open my $fh_pipe, '-|', $cmd;

    my $seq;
    while ( my $line = <$fh_pipe> ) {
        chomp $line;
        if ( $line =~ /^[\w-]+/ ) {
            $seq .= $line;
        }
    }
    close($fh_pipe);

    return $seq;
}

sub decode_header {
    my $header = shift;

    # S288C.chrI(+):27070-29557|species=S288C
    my $head_qr = qr{
                ([\w_]+)?           # name
                [\.]?               # spacer
                ((?:chr)?[\w-]+)    # chr name
                (?:\((.+)\))?       # strand
                [\:]                # spacer
                (\d+)               # chr start
                [\_\-]              # spacer
                (\d+)               # chr end
            }xi;

    tie my %info, "Tie::IxHash";

    $header =~ $head_qr;
    my $name     = $1;
    my $chr_name = $2;

    if ( defined $name ) {
        %info = (
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        );
        if ( !defined $info{chr_strand} ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '1' ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '-1' ) {
            $info{chr_strand} = '-';
        }
    }
    elsif ( defined $chr_name ) {
        $name = $header;
        %info = (
            chr_name   => $2,
            chr_strand => $3,
            chr_start  => $4,
            chr_end    => $5,
        );
        if ( !defined $info{chr_strand} ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '1' ) {
            $info{chr_strand} = '+';
        }
        elsif ( $info{chr_strand} eq '-1' ) {
            $info{chr_strand} = '-';
        }
    }
    else {
        $name = $header;
        %info = (
            chr_name   => 'chrUn',
            chr_strand => '+',
            chr_start  => undef,
            chr_end    => undef,
        );
    }
    $info{name} = $name;

    # additional keys
    if ( $header =~ /\|(.+)/ ) {
        my @parts = grep {defined} split /;/, $1;
        for my $part (@parts) {
            my ( $key, $value ) = split /=/, $part;
            if ( defined $key and defined $value ) {
                $info{$key} = $value;
            }
        }
    }

    return \%info;
}

#----------------------------#
# change fasta sequence names
#----------------------------#
sub change_name_chopped {
    my $seq_of       = shift;
    my $seq_names    = shift;
    my $head_chopped = shift;
    my $tail_chopped = shift;

    my $new_seq_of = {};
    my $new_names  = [];

    for my $n ( @{$seq_names} ) {
        my ( $chr, $set, $strand ) = string_to_set($n);
        my $start = $set->min;
        my $end   = $set->max;

        if ( $strand eq '+' ) {
            if ( $head_chopped->{$n} ) {
                $start = $start + $head_chopped->{$n};
            }
            if ( $tail_chopped->{$n} ) {
                $end = $end - $tail_chopped->{$n};
            }
        }
        else {
            if ( $head_chopped->{$n} ) {
                $end = $end - $head_chopped->{$n};
            }
            if ( $tail_chopped->{$n} ) {
                $start = $start + $tail_chopped->{$n};
            }
        }

        my $new_set = AlignDB::IntSpan->new;
        $new_set->add_pair( minmax( $start, $end ) );
        my $new_name = $chr . "(" . $strand . "):" . $new_set->runlist;

        push @{$new_names}, $new_name;
        $new_seq_of->{$new_name} = $seq_of->{$n};
    }

    return ( $new_seq_of, $new_names );
}

1;
