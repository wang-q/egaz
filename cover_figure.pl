#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Bio::Seq;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

use Path::Tiny;
use Number::Format qw(:subs);
use List::Util qw(max sum0);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/lib";
use MyUtil qw(read_sizes);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

alignDB_graph.pl - Generate graph for chromosome coverage in alignDB

=head1 SYNOPSIS

    perl alignDB_graph.pl [options]
      Options:
        --help      -?          brief help message
        --class         STR     GD or GD::SVG
        --width         INT     width + 40 = figure width
        --goal          STR     coverage on target or query, default is [target]

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $yml_file,
    'size|s=s' => \my $size_file,
    'name|n=s' => \( my $name ),
    'class=s'  => \( my $CLASS = "GD" ),
    'width=i'  => \( my $width = 1000 ),
) or HelpMessage(1);

if (! defined $name) {
    $name = path($yml_file)->basename(".yaml", ".yml");                  
}

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Drawing coverage...");

my $yml = LoadFile($yml_file);

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#

my $length_of        = read_sizes($size_file);
my $length_of_genome = sum0( values %{$length_of} );
my @chrs             = keys %{$length_of};

my $largest = max( values %{$length_of} );
my @chr_infos;

# for each chromosome
for my $chr_name (@chrs) {
    my $chr_set    = AlignDB::IntSpan->new;
    my $chr_length = $length_of->{$chr_name};
    $chr_set->add_pair( 1, $chr_length );
    my $cover_set = AlignDB::IntSpan->new;

    if ( exists $yml->{$chr_name} ) {
        $cover_set->add( $yml->{$chr_name} );
    }

    my $coverage = $cover_set->size / $chr_length * 100;
    $coverage = sprintf "%.2f", $coverage;

    my $info = {
        chr_name   => $chr_name,
        chr_length => $chr_length,
        chr_set    => $chr_set,
        cover_set  => $cover_set,
        coverage   => $coverage,

        #band         => $band,
    };

    push @chr_infos, $info;
}

#----------------------------#
# draw pic
#----------------------------#
print "Draw picture...\n";
my $ftr = 'Bio::Graphics::Feature';
{
    my $largest_chr = $ftr->new(
        -start => 1,
        -end   => $largest,
        -name  => $name,
        -type  => 'alignment'
    );
    my $panel = Bio::Graphics::Panel->new(
        -grid        => 1,
        -gridcolor   => 'lightcyan',
        -segment     => $largest_chr,
        -spacing     => 15,
        -width       => $width,
        -pad_top     => 20,
        -pad_bottom  => 20,
        -pad_left    => 20,
        -pad_right   => 20,
        -key_style   => 'none',
        -image_class => $CLASS,
    );

    {    # arrow
        $panel->add_track(
            $largest_chr,
            -glyph      => 'arrow',
            -double     => 1,
            -fgcolor    => 'red',
            -bump       => 0,
            -height     => 10,
            -arrowstyle => 'regular',
            -tick       => 2,
            -linewidth  => 1,
        );
    }

    for my $chr_info (@chr_infos) {
        my $chr_name   = $chr_info->{chr_name};
        my $chr_set    = $chr_info->{chr_set};
        my $chr_length = $chr_info->{chr_length};
        my $cover_set  = $chr_info->{cover_set};
        my $coverage   = $chr_info->{coverage};

        #my $band         = $chr_info->{band};

        my $chr_segment = $ftr->new(
            -segments => [ $chr_set->spans ],
            -name     => $chr_name,
            -type     => 'alignment'
        );

        {    # text
            my $title
                = "$chr_name: " . format_bytes($chr_length) . " bp" . " | Coverage: $coverage%";
            $panel->add_track(
                $chr_segment,
                -glyph        => 'text_in_box',
                -text         => $title,
                -text_bgcolor => 'lightcyan',
                -height       => 10,
                -bgcolor      => 'yellow',
                -text_pad     => 4,
            );
        }

        {    # alignment
            $panel->add_track(
                $chr_segment,
                -glyph     => 'segments',
                -label     => $chr_name,
                -bump      => 0,
                -height    => 10,
                -font      => 'gdSmallFont',
                -linewidth => 1,
                -bgcolor   => 'lightblue',
                -connector => 'solid',
            );

            my $target_segment = $ftr->new(
                -segments => [ $cover_set->spans ],
                -name     => 'coverage',
                -type     => 'alignment'
            );
            $panel->add_track(
                $target_segment,
                -glyph     => 'segments',
                -label     => 'coverage',
                -bump      => 0,
                -height    => 10,
                -font      => 'gdSmallFont',
                -linewidth => 1,
                -bgcolor   => 'lightgreen',
                -fgcolor   => 'lightgreen',
                -connector => 'solid',
            );
        }
        #
        ## bands
        #if ( defined $band ) {
        #    my $band_seg   = {};
        #    my $band_color = {
        #        acen    => 'green',
        #        gneg    => 'cyan',
        #        gpos100 => 'gray',
        #        gpos75  => 'gray',
        #        gpos50  => 'gray',
        #        gpos25  => 'gray',
        #        gvar    => 'purple',
        #        stalk   => 'green',
        #    };
        #    for (@$band) {
        #        my $source   = $_->stain;
        #        my $band_ftr = Bio::Graphics::Feature->new(
        #            -start => $_->start,
        #            -stop  => $_->end,
        #            -name  => $_->name,
        #            -type  => 'band',
        #
        #            #-source => $_->stain,
        #        );
        #        if ( exists $band_seg->{$source} ) {
        #            push @{ $band_seg->{$source} }, $band_ftr;
        #
        #        }
        #        else {
        #            $band_seg->{$source} = [$band_ftr];
        #        }
        #    }
        #
        #    for my $source ( sort keys %$band_seg ) {
        #        my $band_segment = $ftr->new(
        #            -segments => $band_seg->{$source},
        #            -name     => $source,
        #        );
        #
        #        $panel->add_track(
        #            $band_segment,
        #            -glyph   => 'segments',
        #            -label   => $source,
        #            -fgcolor => $band_color->{$source},
        #            -bgcolor => 'white',
        #            -height  => 10,
        #        );
        #    }
        #}

        # line
        $panel->add_track(
            $largest_chr,
            -glyph     => 'line',
            -bump      => 0,
            -height    => 1,
            -linewidth => 1,
            -bgcolor   => 'turquoise',
            -fgcolor   => 'turquoise',
        );
    }

    my $gd = $panel->gd;
    my $type = ( $CLASS eq 'GD' ) ? 'png' : 'svg';
    open my $pic_fh, '>', "$name.$type";
    binmode $pic_fh;
    print {$pic_fh} $gd->$type;
    close $pic_fh;
}

$stopwatch->end_message;
exit;

__END__
