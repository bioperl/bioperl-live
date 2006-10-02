package Bio::Graphics::Glyph::text_in_box;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_text
{
  return "3'";  
}

sub default_text_pad
{
  return 3;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  
  my $font = $self->option('labelfont') || $self->font;
  
  my $text = defined $self->option('text') ? $self->option('text') : $self->default_text();
  my $text_pad = defined $self->option('text_pad') ? $self->option('text_pad') : $self->default_text_pad();
  
  my $width = $font->width * length $text;
  my $height = $font->height;

  my $midY = ($y2+$y1) / 2;

  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($x1,$midY-$height/2-$text_pad);
  $polygon->addPt($x1+$width+2*$text_pad,$midY-$height/2-$text_pad);
  $polygon->addPt($x1+$width+2*$text_pad,$midY+$height/2+$text_pad);
  $polygon->addPt($x1, $midY+$height/2+$text_pad);

  if (defined (my $bgcolor = $self->option('text_bgcolor')))
  {
    $gd->filledPolygon($polygon,$self->factory->translate_color($bgcolor));
  }

  $gd->polygon($polygon,$fg);
      
  $gd->string($font, $x1+$text_pad, $midY-$height/2, $text, $self->fontcolor);
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::text_in_box - The "text in box" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws the specified text in a rectangular box.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -text       The text to draw in the box    3'

  -text_pad   The number of pixels to offset 3
              the box

  -text_bgcolor                             none
              The background color of the box

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
