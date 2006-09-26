package Bio::Graphics::Glyph::pentagram;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub pad_top
{
  my ($self) = @_;
  
  my $font = $self->option('labelfont') || $self->font;
  
  my $pad = $font->height;
  
  if ($self->option('text'))
  {
    $pad *= 2;
  }
  return $pad;
}

sub default_text
{
  return '';  
}

sub default_text_pad_x
{
  return 0;  
}

sub default_text_pad_y
{
  return 3;  
}

sub default_size
{
  return 20;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;
  
  my $size = defined $self->option('size') ? $self->option('size') : $self->default_size();

  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();

  if ($self->option('inverted') == 1)
  {
    $polygon->addPt($x1,$y2);
    $polygon->addPt($x1+$size/2,$y2-$size/2);
    $polygon->addPt($x1,$y2-$size);
    $polygon->addPt($x1+$size, $y2-$size);
    $polygon->addPt($x1+$size, $y2);      
  }
  else
  {
    $polygon->addPt($x1,$y2);
    $polygon->addPt($x1,$y2-$size);
    $polygon->addPt($x1+$size/2,$y2-$size);
    $polygon->addPt($x1+$size, $y2-$size/2);
    $polygon->addPt($x1+$size/2, $y2);      
  }
  
  $gd->filledPolygon($polygon, $bg);
  $gd->polygon($polygon,$fg);

  my $text = defined $self->option('text') ? $self->option('text') : $self->default_text();

  if ($text)
  {
    my $text_pad_x = defined $self->option('text_pad_x') ? $self->option('text_pad_x') : $self->default_text_pad_x();
    my $text_pad_y = defined $self->option('text_pad_y') ? $self->option('text_pad_y') : $self->default_text_pad_y();
    my $font = $self->option('labelfont') || $self->font;
    $gd->string($font, $x1+$text_pad_x, $y2-$size-$text_pad_y-$font->height, $text, $fg);
    
  }
  
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::pentagram - The "pentagram" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a pentagram with the sharp angle pointing right
or,if the 'inverted' option is set to 1, an "inverted" pentagram
(with the sharp angle pointing inwards, not outwards).
There may be an optional text above the glyph.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -size       Width and height of the       20
              glyph

  -text       Text to show                  none

  -text_pad_x Number of pixels between        0
              the left edge of the glyph
              and the start of text

  -text_pad_x Number of pixels between        3
              the pentagram
              and the text

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
