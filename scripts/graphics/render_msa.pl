#!/usr/bin/perl
use strict;
use Bio::Graphics;
use Bio::AlignIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GuessSeqFormat;
use Getopt::Long;

my ($inputfile,$debug);

GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
           'debug' => \$debug,
          );

my $guessed_format = new Bio::Tools::GuessSeqFormat
    (-file=>"$inputfile"
    )->guess;

my $aio = Bio::AlignIO->new(-file   => $inputfile,
                            -format => $guessed_format) or die "parse failed";

my $aln = $aio->next_aln() or die "no alignment";

my $panel = Bio::Graphics::Panel->new(
    -image_class => 'SVG',
    -length    => $aln->length,
    -width     => 800,
    -pad_left  => 150,
    -pad_right => 10,
);

my $full_length = Bio::SeqFeature::Generic->new(
    -start        => 1,
    -end          => $aln->length,
    -display_name => 'fasta alignment',
);
$panel->add_track($full_length,
                  -glyph   => 'arrow',
                  -tick    => 2,
                  -fgcolor => 'black',
                  -double  => 1,
                  -label   => 1,
                 );
 
my $track = $panel->add_track(
                              -glyph       => 'segments',
                              -label       => 1,
                              -connector   => 'dashed',
                              -bgcolor     => 'green',
                              -label_position   => 'alignment_left',
                              -font2color  => 'red'
                             );
 
for my $seqobj ($aln->each_seq) {
    my $seq = $seqobj->seq;
    my @seqs;
    # get alignment positions for seqs
    my $feature = Bio::SeqFeature::Generic->new(
        -display_name => $seqobj->get_nse,
    );
    while ($seq =~ m{([^-]+)}g) {
        # zero-based coords, must adjust accordingly
        $feature->add_sub_SeqFeature(
        Bio::SeqFeature::Generic->new(-start => pos($seq)-length($1)+1,
                                      -end => pos($seq),
                                      -sequence => $1), 'EXPAND');

    }
    $track->add_feature($feature);
}

print $panel->svg;


