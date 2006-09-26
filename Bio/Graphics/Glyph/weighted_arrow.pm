package Bio::Graphics::Glyph::weighted_arrow;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

use Math::Trig;

sub default_weight_size
{
  return 8;  
}

sub default_length
{
  return 20;  
}

sub default_left_alignment
{
  return 1;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;
  
  my $weight_size = defined $self->option('weight_size') ? $self->option('weight_size') : $self->default_weight_size();

  my $length = defined $self->option('length') ? $self->option('length') : $self->default_length();
  my $left_alignment = defined $self->option('left_alignment') ? $self->option('left_alignment') : $self->default_left_alignment();
  
  my ($w_x1, $l_x1);
  if ($left_alignment == 1)
  {
    $w_x1 = $x1+1;
    $l_x1 = $x1;
  }
  else
  {
    $w_x1 = $x2-1-$weight_size;
    $l_x1 = $x2;
    
  }
  
  my $poly_pkg = $self->polygon_package;
  my $polygon   = $poly_pkg->new();
  $polygon->addPt($w_x1,$y2);
  $polygon->addPt($w_x1,$y2-$weight_size);
  $polygon->addPt($w_x1+$weight_size,$y2-$weight_size);
  $polygon->addPt($w_x1+$weight_size,$y2);
  $gd->filledPolygon($polygon,$bg);

  my $angle = deg2rad(30);
  my $dx = 6;
  my $dy = 4;
  my $midX = $x2-$dx;
  my $midY = $y1+$dy;
  
  $gd->line($l_x1, $y2, $l_x1, $y1+$dy, $fg);

  return unless $left_alignment == 1;
    
  return if ($x2-$x1 <= $weight_size);
  
  $x2 = $x1+$weight_size+1+$length;
  $gd->line($l_x1, $midY, $x2, $midY, $fg);
  $gd->line($x2, $midY, $x2-$dx, $y1, $fg);
  $gd->line($x2, $midY, $x2-$dx, $y1+2*$dy, $fg);
   
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::weighted_arrow - The "weighted arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws an arrow which has is "weighted" by a square
on the left side of the glyph or a "weight" and a vertical line, but
no arrow on its right side. The arrow/line is drawn with the foreground
color, the square - with the background color.

The first mode is the default. To get the second mode, specify
the "left_alignment 0" option.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -weight_size
              Size of the square            8

  -length     Length of the arrow           20

  -left_alignment
              Whether the glyph is drawn    1
              on the left or on the right
              side of the available space.

  -height     Standard option, but          10
              important here

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
