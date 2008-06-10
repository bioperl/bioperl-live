package Bio::Graphics::Glyph::flag;

# $Id$
# Non object-oriented utilities used here-and-there in Bio::Graphics modules

=head1 NAME

Bio::Graphics::Glyph::flag - the "flag" glyph

=cut

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_text
{
  return "ori";  
}

sub default_width
{
  return 20;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;
  
  my $width = $self->option('width') || $self->default_width;
  my $text = $self->option('text') || $self->default_text;
  
  my $oneThirdY = $y1 + ($y2-$y1) / 3;
  my $twoThirdsY = $y1 + 2 * ($y2-$y1) / 3;
 
  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($x1, $y1);
  $polygon->addPt($x1+$width, $oneThirdY);
  $polygon->addPt($x1, $twoThirdsY);

  $gd->polygon($polygon, $fg);
  
  $gd->fillToBorder($x1+$width/2, $oneThirdY, $fg, $bg);  
  
  $gd->line($x1, $y1, $x1, $y2, $fg);
  
  my $font = $self->option('labelfont') || $self->font;

  $gd->string($font, $x1 + 3, $twoThirdsY-3, $text, $self->fontcolor);  
  
}

1;

__END__

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a flag with a text next to it.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -text       Text to draw next to the flag  ori

  -width      Width of the flag               20

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

Vsevolod (Simon) Ilyushchenko E<lt>simonf@cshl.eduE<gt>.

Copyright (c) 2004 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
