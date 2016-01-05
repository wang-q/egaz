package MyUtil;
use strict;
use warnings;
use autodie;

use 5.010001;

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
            string_to_set set_to_string change_strand read_sizes revcom exec_cmd run_sparsemem get_seq_faidx
            get_size_faops decode_header encode_header change_name_chopped sort_cc
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub string_to_set {
    my $node = shift;

    my ( $chr_part, $runlist ) = split /:/, $node;

    my $head_qr = qr{
        (?:(?P<name>[\w_]+)\.)?    
        (?P<chr_name>[\w-]+)        
        (?:\((?P<chr_strand>.+)\))?  
    }xi;
    $chr_part =~ $head_qr;

    my $chr    = $+{chr_name};
    my $strand = "+";
    if ( defined $+{chr_strand} ) {
        $strand = $+{chr_strand};
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

sub run_sparsemem {
    my $file   = shift;
    my $genome = shift;
    my $length = shift || 20;

    my $cmd = sprintf "sparsemem -maxmatch -F -l %d -b -n -k 3 -threads 3 %s %s", $length, $genome,
        $file;
    my $result = `$cmd`;

    return $result;
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
    close $fh_pipe;

    return $seq;
}

sub get_size_faops {
    my $file = shift;

    my $cmd = sprintf "faops size %s", $file;
    my @lines = grep {defined} split /\n/, `$cmd`;

    tie my %length_of, "Tie::IxHash";
    for (@lines) {
        my ( $key, $value ) = split /\t/;
        $length_of{$key} = $value;
    }

    return \%length_of;
}

sub decode_header {
    my $header = shift;

    # S288C.chrI(+):27070-29557|species=S288C
    my $head_qr = qr{
        (?:(?P<name>[\w_]+)\.)?    
        (?P<chr_name>[\w-]+)        
        (?:\((?P<chr_strand>.+)\))? 
        [\:]                        # spacer
        (?P<chr_start>\d+)    
        [\_\-]                      # spacer
        (?P<chr_end>\d+)        
    }xi;

    tie my %info, "Tie::IxHash";

    $header =~ $head_qr;
    my $name     = $1;
    my $chr_name = $2;

    if ( defined $name or defined $chr_name ) {
        %info = (
            name       => $name,
            chr_name   => $chr_name,
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
            name       => undef,
            chr_name   => 'chrUn',
            chr_strand => '+',
            chr_start  => undef,
            chr_end    => undef,
        );
    }

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

sub encode_header {
    my $info           = shift;
    my $only_essential = shift;

    my $header;
    $header .= $info->{name};
    $header .= "." . $info->{chr_name};
    $header .= "(" . $info->{chr_strand} . ")";
    $header .= ":" . $info->{chr_start};
    $header .= "-" . $info->{chr_end};

    # additional keys
    if ( !$only_essential ) {
        my %essential = map { $_ => 1 } qw{name chr_name chr_strand chr_start chr_end seq full_seq};
        my @parts;
        for my $key ( sort keys %{$info} ) {
            if ( !$essential{$key} ) {
                push @parts, $key . "=" . $info->{$key};
            }
        }
        if (@parts) {
            my $additional = join ";", @parts;
            $header .= "|" . $additional;
        }
    }

    return $header;
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


sub sort_cc {
    my @cc = @_;

    # sort by chromosome order within cc
    for my $c (@cc) {
        my @unsorted = @{$c};

        # start point on chromosomes
        @unsorted = map { $_->[0] }
            sort { $a->[1] <=> $b->[1] }
            map { /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @unsorted;

        # chromosome name
        @unsorted = map { $_->[0] }
            sort { $a->[1] cmp $b->[1] }
            map { /([\w.]+)\(.\)\:/; [ $_, $1 ] } @unsorted;

        $c = [@unsorted];
    }

    # sort by first node's chromosome order between cc
    @cc = map { $_->[0] }
        sort { $a->[1] <=> $b->[1] }
        map { $_->[0] =~ /[\w.]+\(.\)\:(\d+)/; [ $_, $1 ] } @cc;

    @cc = map { $_->[0] }
        sort { $a->[1] cmp $b->[1] }
        map { $_->[0] =~ /([\w.]+)\(.\)\:/; [ $_, $1 ] } @cc;

    # sort by nodes number between cc
    @cc = sort { scalar @{$b} <=> scalar @{$a} } @cc;

    return @cc;
}

1;
