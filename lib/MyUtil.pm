package MyUtil;
use strict;
use warnings;
use autodie;

use Carp;
use Path::Tiny;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            read_sizes read_fasta revcom exec_cmd
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub read_sizes {
    my $file       = shift;
    my $remove_chr = shift;

    my @lines = path($file)->lines( { chomp => 1 } );
    my %length_of;
    for (@lines) {
        my ( $key, $value ) = split /\t/;
        $key =~ s/chr0?// if $remove_chr;
        $length_of{$key} = $value;
    }

    return \%length_of;
}

sub read_fasta {
    my $filename = shift;

    my @lines = path($filename)->lines( { chomp => 1,});

    my @seq_names;
    my %seqs;
    for my $line (@lines) {
        if ( $line =~ /^\>[\w:-]+/ ) {
            $line =~ s/\>//;
            push @seq_names, $line;
            $seqs{$line} = '';
        }
        elsif ( $line =~ /^[\w-]+/ ) {
            $line =~ s/[^\w-]//g;
            my $seq_name = $seq_names[-1];
            $seqs{$seq_name} .= $line;
        }
        else {    # Blank line, do nothing
        }
    }

    return ( \%seqs, \@seq_names );
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
