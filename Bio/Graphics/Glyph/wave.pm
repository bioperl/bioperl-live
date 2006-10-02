package Bio::Graphics::Glyph::wave;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_spread
{
  return 0.3;  
}

sub default_radius
{
  return 5;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $spread = defined $self->option('spread') ? $self->option('spread') : $self->default_spread();
  
  my $fg = $self->fgcolor;

  my $height = ($y2-$y1)/2;
  my $midY = $y1 + $height;
  
  if ($self->option('circle') == 1)
  {
    my $radius = defined $self->option('radius') ? $self->option('radius') : $self->default_radius();
    $gd->ellipse($x1+$radius,$midY,2*$radius,2*$radius,$fg);
    $x1 = $x1+2*$radius;
  }
  
  if ($self->option('line') == 1)
  {
    if ($x1 < $x2)
    {
      $gd->line($x1,$midY,$x2,$midY,$fg);
    }
    return;
  }
    
  my ($oldX, $oldY);
  foreach my $x ($x1..$x2)
  {
    my $y = -$height * sin ($spread * ($x-$x1))+$midY;
    if ($x>$x1)
    {
      $gd->line($oldX,$oldY,$x,$y,$fg);
    }
    $oldX=$x;
    $oldY=$y;
  }
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::wave - The "wave" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a sine wave with an optional circle in the beginning.
The wave can also become a straight line.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -spread     The "spread" of the sine curve 0.3
              Values from 0.1 to 0.5 look best

  -line       Whether  to draw a line         0
              instead of a wave (1 or 0)

  -circle     Whether to draw a circle        0
              in the left corner (1 or 0)

  -radius     The radius of the circle        5
              if present

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
