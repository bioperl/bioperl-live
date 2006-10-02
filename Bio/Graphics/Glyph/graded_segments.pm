package Bio::Graphics::Glyph::graded_segments;
#$Id$

use strict;
use base qw(Bio::Graphics::Glyph::minmax Bio::Graphics::Glyph::merge_parts);

# override draw method to calculate the min and max values for the components
sub draw {
  my $self = shift;

  # bail out if this isn't the right kind of feature
  # handle both das-style and Bio::SeqFeatureI style,
  # which use different names for subparts.
  my @parts = $self->parts;
  @parts    = $self if !@parts && $self->level == 0;
  return $self->SUPER::draw(@_) unless @parts;

  my ($min_score,$max_score) = $self->minmax(\@parts);

  return $self->SUPER::draw(@_)
    unless defined($max_score) && defined($min_score)
      && $min_score < $max_score;

  my $span = $max_score - $min_score;

  # allocate colors
  my $fill   = $self->bgcolor;
  my ($red,$green,$blue) = $self->panel->rgb($fill);

  @parts = $self->merge_parts(@parts) if $self->option('merge_parts');

  foreach my $part (@parts) {
    my $s = eval { $part->feature->score };
    unless (defined $s) {
      $part->{partcolor} = $fill;
      next;
    }
    my ($r,$g,$b) = $self->calculate_color($s,[$red,$green,$blue],$min_score,$span);
    my $idx      = $self->panel->translate_color($r,$g,$b);
    $part->{partcolor} = $idx;
  }
  $self->SUPER::draw(@_);
}

sub calculate_color {
  my $self = shift;
  my ($s,$rgb,$min_score,$span) = @_;
  return map { 255 - (255-$_) * min(max( ($s-$min_score)/$span, 0), 1) } @$rgb;
}

sub min { $_[0] < $_[1] ? $_[0] : $_[1] }
sub max { $_[0] > $_[1] ? $_[0] : $_[1] }

sub subseq {
  my $class = shift;
  my $feature = shift;
  return $feature->segments        if $feature->can('segments');
  return $feature->sub_SeqFeature  if $feature->can('sub_SeqFeature');
  return;
}

# synthesize a key glyph
sub keyglyph {
  my $self = shift;

  my $scale = 1/$self->scale;  # base pairs/pixel

  # two segments, at pixels 0->50, 60->80
  my $offset = $self->panel->offset;

  my $feature =
    Bio::Graphics::Feature->new(
				-segments=>[ [ 0*$scale +$offset,20*$scale+$offset],
					     [ 30*$scale +$offset,50*$scale+$offset],
					     [60*$scale+$offset, 80*$scale+$offset]
					   ],
				-name => $self->option('key'),
				-strand => '+1');
  ($feature->segments)[0]->score(10);
  ($feature->segments)[1]->score(50);
  ($feature->segments)[2]->score(100);
  my $factory = $self->factory->clone;
  $factory->set_option(label => 1);
  $factory->set_option(bump  => 0);
  $factory->set_option(connector  => 'solid');
  return $factory->make_glyph($feature);
}

# component draws a shaded box
sub bgcolor { 
  my $self = shift;
  return defined $self->{partcolor} ? $self->{partcolor} : $self->SUPER::bgcolor;
}
sub fgcolor {
  my $self = shift;
  return $self->SUPER::fgcolor unless $self->option('vary_fg');
  return defined $self->{partcolor} ? $self->{partcolor} : $self->SUPER::fgcolor;
}

1;

=head1 NAME

Bio::Graphics::Glyph::graded_segments - The "graded_segments" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is identical to the "alignment" glyph, and is used for
drawing features that consist of discontinuous segments.  The
color intensity of each segment is proportionate to the score.

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

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option      Description                   Default
  ------      -----------                   -------

  -max_score  Maximum value of the	    Calculated
              feature's "score" attribute

  -min_score  Minimum value of the          Calculated
              feature's "score" attribute

  -vary_fg    Vary the foreground color as  0 (false)
              well as the background

  -merge_parts                             0 (false)
              Whether to simplify the
              alignment at low magnification

  -max_gap    Do not merge across gaps     Calculated
              that exceed this threshold


If max_score and min_score are not specified, then the glyph will
calculate the local maximum and minimum scores at run time.  Since
many scoring functions are exponential you may wish to take the log of
your scores before passing them to this glyph.


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


=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::cds>,
L<Bio::Graphics::Glyph::crossbox>,
L<Bio::Graphics::Glyph::diamond>,
L<Bio::Graphics::Glyph::dna>,
L<Bio::Graphics::Glyph::dot>,
L<Bio::Graphics::Glyph::ellipse>,
L<Bio::Graphics::Glyph::extending_arrow>,
L<Bio::Graphics::Glyph::generic>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::heterogeneous_segments>,
L<Bio::Graphics::Glyph::line>,
L<Bio::Graphics::Glyph::pinsertion>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::rndrect>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::ruler_arrow>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::translation>,
L<Bio::Graphics::Glyph::triangle>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
