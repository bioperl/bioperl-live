package Bio::Graphics::Glyph::ex;

use strict;
use base 'Bio::Graphics::Glyph::generic';

# override draw_component to draw a crossed box rather than empty
sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  my ($left,$top) = @_;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_);

  #if widthless
  if($self->option('point')){
    my $arm = int($self->height/2);
    my $minx    = $x2 > $x1 ? $x1 : $x2;
    my $centerx = abs($x2 - $x1) + $minx;
    my $miny    = $y2 > $y1 ? $y1 : $y2;
    my $centery = abs($y2 - $y1) + $miny;
    $gd->line($centerx-$arm, $centery-$arm, $centerx+$arm, $centery+$arm, $fg);
    $gd->line($centerx-$arm, $centery+$arm, $centerx+$arm, $centery-$arm, $fg);
    return;
  } else {
    $gd->line($x1,$y1,$x2,$y2,$fg);
    $gd->line($x1,$y2,$x2,$y1,$fg);
  }
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::crossbox - The "crossbox" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This is a box with an 'X' inside glyph.

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

Allen Day E<lt>day@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
