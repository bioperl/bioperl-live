package Bio::Graphics::Glyph::ragged_ends;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  my $bg = $self->option('bgcolor');
  my ($left,$top) = @_;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $zig = $self->option('zig') || 4;
  my $zag = $self->option('zag') || 4;

  my $polygon = GD::Polygon->new;

  $polygon->addPt($x1, $y1);
  my $yoff = $y1 + $zig;
  my $i = 1;
  if ($self->option("ragged_left")) {
      while ($yoff <= $y2) {
	  $polygon->addPt( $x1 + ($i * $zag),
			   $yoff );
	  $i = !$i;
	  $yoff += $zig;
      }
  }
  $polygon->addPt($x1, $y2);

  $polygon->addPt($x2, $y2);
  $yoff = $y2 - $zig;
  $i = 1;
  if ($self->option("ragged_right")) {
      while ($yoff >= $y1) {
	  $polygon->addPt( $x2 - ($i * $zag),
			   $yoff );
	  $i = !$i;
	  $yoff -= $zig;
      }
  }
  $polygon->addPt($x2, $y1);

  $polygon->addPt($x1, $y1); # close the polygon

  $gd->polygon($polygon, $fg);
  $gd->filledPolygon($polygon, $bg) if $bg;
}

1;

=head1 NAME

Bio::Graphics::Glyph::ragged_ends - The "ragged ends" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is identical to the "box" glyph except that it will draw the
subparts of features that contain subfeatures.  The subparts are not
connected -- use the "segments" glyph for that.  "Generic" is the
default glyph used when not otherwise specified.

=head2 OPTIONS

This glyph provides two extra options to control whether the right
and/or left ends of the drawn box are to be drawn "raggedly" with
zigzags instead of vertical lines.

  Option        Values    Default
  -raggedleft    0 | 1       1
  -raggedright   0 | 1       1
  -zig           > 3         4
  -zag           > 3         4

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

  -pad_top      Top padding                    0

  -pad_bottom   Bottom padding                 0

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

  -hilite       Highlight color                undef (no color)

-pad_top and -pad_bottom allow you to insert some blank space between
the glyph's boundary and its contents.  This is useful if you are
changing the glyph's height dynamically based on its feature's score.

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
L<Bio::Graphics::Glyph::xyplot>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Aaron J Mackey E<lt>amackey@pcbi.upenn.eduE<gt>.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
