package Bio::Graphics::Glyph::diamond;
# DAS-compatible package to use for drawing a colored diamond

use strict;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;

  # find the center and vertices
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  my $xmid = ($x1+$x2)/2;
  my $ymid = ($y1+$y2)/2;

  my $h = $self->option('height')/2;
  $x1 = $xmid - $h;
  $x2 = $xmid + $h;
  $y1 = $ymid - $h;
  $y2 = $ymid + $h;

  # now draw the diamond
  $gd->line($x1,$ymid,$xmid,$y1,$fg);
  $gd->line($xmid,$y1,$x2,$ymid,$fg);
  $gd->line($x2,$ymid,$xmid,$y2,$fg);
  $gd->line($xmid,$y2,$x1,$ymid,$fg);

  if (my $c = $self->bgcolor) {
    $gd->fillToBorder($xmid,$ymid,$fg,$c);
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::diamond - The "diamond" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a diamond of fixed size, positioned in the center of
the feature.  The height and width of the diamond are set by the
"height" option.

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

=head1 BUGS

If the feature is wider than a point, then the label and description
are placed where the feature's boundary is, and not where the diamond
is drawn.

=head1 SEE ALSO

L<Ace::Sequence>, L<Ace::Sequence::Feature>, L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>, L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,

=head1 AUTHOR

Lincoln Stein <lstein@cshl.org>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
