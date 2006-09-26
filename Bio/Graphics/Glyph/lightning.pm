package Bio::Graphics::Glyph::lightning;

# A lightning bolt glyph to add some pizazz to your displays. Yeow!

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub draw_component {
  my $self = shift;
  my $gd   = shift;
  my $fg   = $self->fgcolor;

  # find the center and vertices
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $h = $self->option('height');
  my $w = $h*0.6;

  my $poly_pkg = $self->polygon_package;
  my $polygon   = $poly_pkg->new();

  # lightning bolt points up or down
  if ($self->option('orient') eq 'N') {
    $y1 = $y1 + $h/4;
    $polygon->addPt($x2,$y1+($h/2));
    $polygon->addPt($x2-($w*0.7),$y1+($h/2));
    $polygon->addPt($x2-($w*0.2),$y1+($h*0.15));
    $polygon->addPt($x2-($w*0.6),$y1+($h*0.15));
    $polygon->addPt($x2,$y1-$h/2);
    $polygon->addPt($x2-($w*0.1),$y1-($h*0.05));
    $polygon->addPt($x2+($w*0.5),$y1-($h*0.05));
    $polygon->addPt($x2,$y1+($h/2));
  }
  else {
    $y1 = $y1 + $h/2;
    $polygon->addPt($x1,$y1-($h/2));
    $polygon->addPt($x1+($w*0.7),$y1-($h/2));
    $polygon->addPt($x1+($w*0.2),$y1-($h*0.15));
    $polygon->addPt($x1+($w*0.6),$y1-($h*0.15));
    $polygon->addPt($x1,$y1+$h/2);
    $polygon->addPt($x1+($w*0.1),$y1+($h*0.05));
    $polygon->addPt($x1-($w*0.5),$y1+($h*0.05));
    $polygon->addPt($x1,$y1-($h/2));
}

  # Have to draw TWO polygons for fills in order to get an outline
  # because filledPolygon in GD croaks with extra parameters (and
  # doesn't support drawing of stroke anyways).
  if (my $c = $self->bgcolor) {
      $gd->filledPolygon($polygon,$c);
      $gd->polygon($polygon,$fg);
  } else {
    $gd->polygon($polygon,$fg);
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::lightning - The "lightning" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a lightning bolt of specified height with relative
width, with the point of the lightning bolt centered on the
feature. The height of the bolt is specified by the "height"
option. Due to the complexity of this glyph, it doesn't resolve well
with heights less than 11 pixels.

This glyph was designed to indicate point mutations on a nucleotide or
protein backbone.

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



The following options are specific to this Glyph.

  Option      Description                      Default
  ------      -----------                      -------  
  -orient     direction of lightning bolt      N


=head1 BUGS

No reported bugs.

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

Todd Harris E<lt>harris@cshl.orgE<gt>

Copyright (c) 2004 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
