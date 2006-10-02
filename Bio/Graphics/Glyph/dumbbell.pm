package Bio::Graphics::Glyph::dumbbell;
# DAS-compatible package to use for drawing a line of repeating shapes

use strict;
use base qw(Bio::Graphics::Glyph::generic);

use Math::Trig;

sub default_shape_size
{
  return 10;  
}

sub default_shape
{
  return 'square';  
}

sub draw_end_shape
{
  my ($self, @args) = @_;
  my $shape = $self->option('end_shape') || $self->default_shape();
  my $method = "draw_end_$shape";
  if ($self->can($method))
  {
    return $self->$method(@args);
  }
  else
  {
    return $self->draw_end_square(@args);  
  }
}

sub draw_end_square
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;

  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;

  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($x1,$y1);
  $polygon->addPt($x2,$y1);
  $polygon->addPt($x2,$y2);
  $polygon->addPt($x1, $y2);

  $gd->filledPolygon($polygon,$fg);
  
  return ($x1, $x2);
}

sub draw_end_diamond
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;

  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;

  my $poly_pkg = $self->polygon_package;

  my $midX = ($x1+$x2)/2;
  my $midY = ($y1+$y2)/2;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($x1,$midY);
  $polygon->addPt($midX,$y1);
  $polygon->addPt($x2,$midY);
  $polygon->addPt($midX,$y2);

  $gd->filledPolygon($polygon,$fg);
  
  return ($x1, $x2);  
}

sub translated_polygon
{
  my ($self, $midX, $midY, $scale_factor, @coords) = @_;

  my $poly_pkg = $self->polygon_package;

  my $polygon   = $poly_pkg->new();
  for (my $i=0; $i<(scalar @coords) / 2; $i++)
  {
    $polygon->addPt($coords[2*$i], $coords[2*$i+1]);
  }
  
  $polygon->scale($scale_factor, $scale_factor);
  $polygon->offset($midX, $midY);
	
	return $polygon;
}

sub draw_end_star
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;
  
  #my @coords = (95, -31, -58, 81, 0, -100, 58, 81, -95, -31);
  my @coords1 = (31, 42, 31, -42, -49, -30, -38, 0, -49, 30);
  my @coords2 = (100, 0, -81, 59, 31, -95, 31, 95, -81, -58);

  my $star_size = 190;
  
  my $scale_factor = $shape_size / $star_size;
  
  my ($midX, $midY) = ($x1+$shape_size/2, $y1+$shape_size/2);
 
	$gd->filledPolygon($self->translated_polygon($midX, $midY, $scale_factor,  @coords1), $fg);
	$gd->filledPolygon($self->translated_polygon($midX, $midY, $scale_factor, @coords2), $fg);

  return ($midX, $midX);
}
 
sub draw_end_tree
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;

  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;
  
  my $trunk_width = $shape_size/4;

  my $midX = ($x1+$x2)/2;
  my $midY = ($y1+$y2)/2;

  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($midX-$trunk_width/2,$midY);
  $polygon->addPt($midX+$trunk_width/2,$midY);
  $polygon->addPt($midX+$trunk_width/2,$y2);
  $polygon->addPt($midX-$trunk_width/2,$y2);

  $gd->filledPolygon($polygon, $fg);
  
  $self->filled_oval($gd, $x1, $y1, $x2, $y1+2*$shape_size/3, $fg, $fg);
  
  return ($midX, $midX);
}

sub draw_end_clover
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;

  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;
  
  my $trunk_width = $shape_size/4;

  my $midX = ($x1+$x2)/2;

  my $poly_pkg = $self->polygon_package;
  
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($midX-$trunk_width/2,$y1+0.4*$shape_size);
  $polygon->addPt($midX+$trunk_width/2,$y1+0.4*$shape_size);
  $polygon->addPt($midX+$trunk_width/2,$y2);
  $polygon->addPt($midX-$trunk_width/2,$y2);

  $gd->filledPolygon($polygon, $fg);
  
  my $radius = $shape_size / 4.3;
  
  $self->filled_oval($gd, $midX-$radius, $y1, $midX+$radius, $y1+2*$radius, $fg, $fg);
  $self->filled_oval($gd, $x1, $y1+1.3*$radius, $x1+2*$radius, $y1+3.3*$radius, $fg, $fg);
  $self->filled_oval($gd, $x2-2*$radius, $y1+1.3*$radius, $x2, $y1+3.3*$radius, $fg, $fg);
  
  return ($midX, $midX);
}

sub draw_end_bubble
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;
  
  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;

  my $midX = ($x1+$x2)/2;
  my $midY = ($y1+$y2)/2;
  
  my $bubble_text = defined $self->option('bubble_text') ? $self->option('bubble_text') : "Text";

  my $font = $self->option('labelfont') || $self->font;
  my $bubble_text_length = $font->width * length($bubble_text);
  my $bubble_text_x = $midX -  $bubble_text_length / 2;
  my $bubble_text_y = $midY - $font->height / 2;
  
  $gd->string($font, $bubble_text_x, $bubble_text_y, $bubble_text, $self->fontcolor);
  
  my $oval_width = $bubble_text_length * 1.414;
  my $oval_height = $font->height * 1.414;

  $self->oval($gd, $midX-$oval_width/2, $midY-$oval_height/2, $midX+$oval_width/2, $midY+$oval_height/2);

  return ($midX-$oval_width/2, $midX+$oval_width/2);
}

sub draw_end_arrow
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg, $antiparallel) = @_;
  
  my $x2 = $x1 + $shape_size;
  my $y2 = $y1 + $shape_size;
  
  my $angle = deg2rad(30);
  my $dx = 2*$shape_size*cos($angle)/5;
  my $dy = 2*$shape_size*sin($angle)/5;
  my $midX = $x2-$dx;
  my $midY = ($y1+$y2)/2;

  $gd->line($x1, $midY, $x2, $midY, $fg);
  if ($antiparallel)
  {
    $gd->line($x1, $midY, $x1+$dx, $midY-$dy, $fg);
    $gd->line($x1, $midY, $x1+$dx, $midY+$dy, $fg);
  }  
  else
  {
    $gd->line($x2, $midY, $x2-$dx, $midY-$dy, $fg);
    $gd->line($x2, $midY, $x2-$dx, $midY+$dy, $fg);
  }  
  return ($x1, $x2);
}

sub draw_end_wave
{
  my ($self, $gd, $x1, $y1, $shape_size, $fg) = @_;
  
  my $x2 = $x1 + $shape_size;
  
  #Make the heigh constant
  my $y2 = $y1 + $shape_size/2;
  $y1 = $y2-10;
  
  my $step = $shape_size/6;
  $gd->line($x1, $y2, $x1+$step, $y1, $fg);
  $gd->line($x1+$step, $y1, $x1+2*$step, $y2, $fg);
  $gd->line($x1+2*$step, $y2, $x1+3*$step, $y1, $fg);
  $gd->line($x1+3*$step, $y1, $x1+4*$step, $y2, $fg);
  $gd->line($x1+4*$step, $y2, $x1+5*$step, $y1, $fg);
  $gd->line($x1+5*$step, $y1, $x1+6*$step, $y2, $fg);
  return ($x1, $x2);
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  
  my $shape_size = defined $self->option('shape_size') ? $self->option('shape_size') : $self->default_shape_size;

  # find the center and vertices
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  if ($x2-$x1 < $shape_size)
  {
    return $self->SUPER::draw_component($gd, @_);  
  }
  
  my $midX = ($x2-$x1) / 2 + $x1;
  my $midY = ($y2-$y1) / 2 + $y1;
  my $startY = $midY - $shape_size/2;
  
  my $antiparallel = $self->option('antiparallel');
  
  #We need to store the bounds of the shapes drawn because the connecting line will have
  #different length depending on them.
  my ($leftX1, $leftX2) = $self->draw_end_shape($gd, $x1, $startY, $shape_size, $fg);
  my ($rightX1, $rightX2) = $self->draw_end_shape($gd, $x2-$shape_size, $startY, $shape_size, $fg, $antiparallel);

  if ($self->option('arc') == 1)
  {
    #Draw an arc of an ellipse relative to the midpoint between shapes
    #whose center is at (0, -q) and which intersects with the X axis at (p,0) and (-p, 0).
    my $p = ($rightX1 - $leftX2) / 2;
    my $q = $shape_size/2;
    
    my $c = 2 * $p / sqrt(3);
    my $d = 2 * $q;
    my $b = $q - $d;
    my $angle = atan2(sqrt(3), 1);
    my $deg = rad2deg($angle);
    $gd->arc($leftX2+$p,$midY+$q,2*$c,2*$d,270-$deg,270+$deg,$self->factory->translate_color('black'));
  }  
  else
  {
    $gd->line($leftX2,$midY,$rightX1,$midY,$fg);  
  }
  
  if (my $caption = $self->option('caption'))
  {
    my $font = $self->option('labelfont') || $self->font;
    my $midX = ($x2-$x1-2*$shape_size)/2+$x1+$shape_size;
    my $startCaption = $midX - $font->width * length($caption) / 2;
    $gd->string($font, $startCaption, $midY-$font->height, $caption, $self->fontcolor);
  }
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::dumbbell - A glyph that draws a "dumbbell" with the same shapes on both ends. 

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a "dumbbell" with the same shapes on both ends. 

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                       Default
  ------      -----------                       -------

  -shape_size The size of the shape               10
              on both ends.

  -end_shape  One of 'square', 'diamond',         square
              'tree', 'clover', 'star',
              'bubble', 'arrow', 'wave'

  -bubble_text The text to show in the bubble     Text
                if the bubble option is chosen
                above (shape_size is then ignored)

  -antiparallel Whether the right arrow               0
                is reversed

  -arc        Whether the shapes are               0
              connected by an arc
              (a straight line is the default).

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
