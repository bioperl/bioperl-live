package Bio::Graphics::Glyph::merged_alignment;

# $Id$

# this glyph acts like graded_segments but the bgcolor of each segment is
# more configurable.  Supply a list of colors and corresponding
# bins.  Each bin is a range of x-y scores, where score n is > x and <= y
# e.g.
# [ALIGNMENT]
# feature    = alignment
# bins       = 0-50 50-70 70-90 90-100 
# bincolors  = white powderblue cornflowerblue blue  

# for sematic zooming (at lower magnification), to
# reduce visual complexity of alignment
# [ALIGNMENT:20000] # segment length >= 20000 
# feature     = alignment
# merge_parts = 1
# max_gap     = 500 # do not merge across gaps > 500 bp 


use strict;
use base qw(Bio::Graphics::Glyph::graded_segments);

use constant COLORS => "lightgrey powderblue cornflowerblue blue";

# override draw method
sub draw {
  my $self = shift;

  # bail out if this isn't the right kind of feature
  # handle both das-style and Bio::SeqFeatureI style,
  # which use different names for subparts.
  my @parts = $self->parts;
  @parts    = $self if !@parts && $self->level == 0;
  return $self->SUPER::draw(@_) unless @parts;

  my $cols = $self->option('bincolors') || COLORS;
  my @cols = split /\s+/, $cols;
  my $bins = $self->option('bins');
  my @bins = $bins ? split /\s+/, $bins : $self->get_bins(\@parts, @cols);
  my %color;
  @color{@bins} = @cols;

  @parts = $self->merge_parts(@parts) if $self->option('merge_parts');

  # figure out the colors
  for my $part (@parts) {
      my ($bin) = grep { $part->in_range($_) } @bins;
    
      my $idx   = $bin ? $self->panel->translate_color($color{$bin}) 
	  : $self->panel->translate_color('white');
      $part->{partcolor} = $idx;
  }
  
  $self->{parts} = \@parts;

  $self->SUPER::draw(@_);
}

sub in_range {
    my $self = shift;
    my $range = shift;
    my ($low,$high) = split '-', $range;
    my $s = $self->score || shift;
    return 1 if $s > $low && $s <= $high;
    return 0;    
}

# overide background method to paint glyph white as
# a last resort
sub bgcolor {
    my $self = shift;
    return $self->{partcolor} || 'white';
}

# used if bins are not defined in the configuration
# makes equal sized bins corresponding to the number of
# colors specified
sub get_bins {
    my $self  = shift;
    my $parts = shift;
    my $cols  = @_;
    my ($min,$max) = $self->minmax($parts);
    my $range = $max - $min;
    return ($max) if $range == 0;
    my $increment = $range/$cols;
    
    my ($score,@bins) = $min;
    until ($score >= $max) {
	my $range = "$score-";
	$score += $increment;
	$range .= $score;
	push @bins, $range;
    }
    
    return @bins;
}


# synthesize a key glyph
sub keyglyph {
  my $self = shift;

  my $scale = 1/$self->scale;  # base pairs/pixel

  # two segments, at pixels 0->50, 60->80
  my $offset = $self->panel->offset;
  
  my $feature =
    Bio::Graphics::Feature->new(
				-segments=>[ [ 0*$scale +$offset,25*$scale+$offset],
					     [ 25*$scale +$offset,50*$scale+$offset],
					     [ 50*$scale+$offset, 75*$scale+$offset]
					   ],
				-name => $self->option('key'),
				-strand => '+1');

  my @scores = $self->example_scores;
  my @segments = $feature->segments;
  for ($feature->segments) {
      $_->score(shift @scores);
  }

  my $factory = $self->factory->clone;
  $factory->set_option(label => 1);
  $factory->set_option(bump  => 0);
  $factory->set_option(connector  => 'solid');
  my $glyph = $factory->make_glyph(0,$feature);
}

sub example_scores {
    my $self = shift;
    my $bins = $self->option('bins');

    if ($bins) {
	my @bins = split /\s+/, $bins;
	$bins[0] =~ s/(\S+)\-\S+/$1/;
	$bins[-1] =~ s/\S+\-(\S+)/$1/;
	my $mid  = $bins[0] + ($bins[-1] - $bins[0])/2;
    
	return ($bins[0], $mid, $bins[-1]);
    }
    if ($self->option('min_score') || $self->option('max_score')) {
	my ($min,$max) = $self->minmax;
	my $mid  = $min + ($max - $min)/2;
	return($min,$mid,$max);
    }
    
    return (0,50,100);
}

1;

=pod

=head1 NAME

Bio::Graphics::Glyph::merged_alignment - The "merged_alignment" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph acts like graded_segments but the bgcolor of segments 
(sub-feature) is controlled by binned scores.  It also supports
semantic zooming to optimize glyph drawing for larger sequence
displays.

=head2 OPTIONS

The following options are standard among all Glyphs.  See
L<Bio::Graphics::Glyph> for a full explanation.

  Option      Description                      Default
  ------      -----------                      -------

  -fgcolor      Foreground color	       black

  -outlinecolor	Synonym for -fgcolor

  -bgcolor      Background color               turquoise

  -fillcolor    Synonym for -bgcolor

  -linewidth    Line width                     1

  -height       Height of glyph		       10

  -font         Glyph font		       gdSmallFont

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

  -hilite       Highlight color                undef (no color)


In addition, the merged-alignment glyph recognizes the following
glyph-specific options:

  Option      Description                  Default
  ------      -----------                  -------

  -max_score  Maximum value of the         Calculated
              feature's "score" attribute

  -min_score  Minimum value of the         Calculated
              feature's "score" attribute

  -bincolors  Colors assigned to bins      lightgrey powderblue cornflowerblue blue
              (in order)

  -bins       Bins to which scores are     Calculated
              assigned

  -merge_parts                             0 (false)
              Whether to simplify the 
              alignment at low magnification

  -max_gap    Do not merge across gaps     Calculated
              that exceed this threshold


If max_score and min_score are not specified, then the glyph will
calculate the local maximum and minimum scores at run time.

If the bins are not specified, they will be calculated
based on the number of colors assigned and the local
(or user-specified) minimum and maximum scores.
Calculated bins are equal in size.  

User-specified bins are expressed as ranges,

  bins  = 0-50 50-70 70-90 90-100

where each range means greater than the lower number and
less than or equal to the higher number.


=head2 Simplifying the display of alignment features for large segments

The "merge_parts" option is used for semantic zooming.
Specifically, if features are small and dense, they
will not be displayed very well for large segments and the 
color-coding will be lost.  If merge-parts is set to a
true value, adjacent alignment parts will be merged until a gap
exceeding a calculated or user-specified value is encountered. 
Unless specified, the maximum gap allowed for merging adjacent features is
calculated as (L/10000)*(L/500), where L = the length of the sequence
displayed in the browser.  The exponentially increasing gap threshold
allows more aggressive merging of alignment features as the size of
the displayed sequence grows larger.

The score of the merged feature is calculated as a weighted average.
For example, consider two adjacent HSPs that are each 400 bp in 
length and have scores of 60% and 70%.  If the merge_parts option
is set to a true value, the two HSPs would be merged in the display to
a single 800 bp alignment block with an average score of 65%.

The merge_parts option is turned off by default.

=head2 SAMPLE CONFIGURATION

Sample gbrowse configuration stanzas for an alignment feature
using this glyph.  The scores are assumed to be expressed 
as percent identity (0-100).

 # base configuration
 [BLASTZ]
 feature      = blastz_alignment
 glyph        = merged_alignment
 bincolors    = #A0A0A0 powderblue cornflowerblue blue
 bins         = 60-70 70-80 80-90 90-100
 category     = Sequence Similarity Tracks
 height       = 6
 bump         = 1
 label        = 1
 fgcolor      = black
 key          = BLASTZ

Semantic zooming with defined maximum gap between
merged features for different zoom levels

 # if the displayed segment is >= 20000 in length,
 # use the merge_parts option to simplify the alignment
 # display
 [BLASTZ:20000]
 feature      = blastz_alignment
 merge_parts  = 1
 max_gap      = 50 # do not merge across gaps > 50 bp

 # if the displayed segment is >= 50000 in length
 [BLASTZ:50000]
 feature      = blastz_alignment
 merge_parts  = 1
 max_gap      = 500 # do not merge across gaps > 500 bp

--OR-- 

Semantic zooming with dynamically calculated maximum
gap

 # if the displayed segment is >= 20000 in length,
 [BLASTZ:20000]
 feature      = blastz_alignment
 merge_parts  = 1

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Sheldon McKay E<lt>mckays@cshl.eduE<gt>

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
