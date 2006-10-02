package Bio::Graphics::Glyph::christmas_arrow;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

use Math::Trig;

sub default_radius
{
  return 4;  
}

sub default_length
{
  return 20;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  
  my $radius = defined $self->option('radius') ? $self->option('radius') : $self->default_radius ();
  
  $gd->filledEllipse($x1+$radius, $y2-$radius, 2*$radius, 2*$radius, $fg);
  
  my $length = defined $self->option('length') ? $self->option('length') : $self->default_length();

  my $angle = deg2rad(30);
  my $dx = 6;
  my $dy = 4;
  my $midX = $x2-$dx;
  my $midY = $y1+$dy;

  $gd->line($x1+$radius, $y2-$radius, $x1+$radius, $y1+$dy, $fg);

  return if ($x2-$x1 <= $radius);
  
  $x2 = $x1+$radius+$length;
  $gd->line($x1+$radius, $midY, $x2, $midY, $fg);
  $gd->line($x2, $midY, $x2-$dx, $y1, $fg);
  $gd->line($x2, $midY, $x2-$dx, $y1+2*$dy, $fg);
   
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::christmas_arrow - The "christmas arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws an arrow which has a circle ("christmas ball")
dangling at one end.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -radius     Radius of the circle          4
              glyph

  -length     Length of the arrow           20

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
