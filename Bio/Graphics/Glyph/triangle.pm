package Bio::Graphics::Glyph::triangle;
# DAS-compatible package to use for drawing a triangle

use strict;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';
use Bio::Graphics::Glyph::generic;

sub pad_left {
  my $self = shift;
  my $left = $self->SUPER::pad_left;
  return $left unless $self->option('point');
  my $extra = $self->option('height')/3;
  return $extra > $left ? $extra : $left;
}

sub pad_right {
  my $self = shift;
  my $right = $self->SUPER::pad_right;
  return $right unless $self->option('point');
  my $extra = $self->option('height')/3;
  return $extra > $right ? $extra : $right;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  my $orient = $self->option('orient') || 'S';

  # find the center and vertices
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  my $xmid = ($x1+$x2)/2;
  my $ymid = ($y1+$y2)/2;

  my ($vx1,$vy1,$vx2,$vy2,$vx3,$vy3);

  #make an equilateral
  my ($p,$q) = ($self->option('height'),($x2-$x1)/2);
  if ($self->option('point')){
    $q = $p/sqrt(3); #2;
    $x1 = $xmid - $q; $x2 = $xmid + $q;
    $y1 = $ymid - $q; $y2 = $ymid + $q;
  }

  if   ($orient eq 'S'){$vx1=$x1;$vy1=$y1;$vx2=$x2;$vy2=$y1;$vx3=$xmid;$vy3=$y2;}
  elsif($orient eq 'N'){$vx1=$x1;$vy1=$y2;$vx2=$x2;$vy2=$y2;$vx3=$xmid;$vy3=$y1;}
  elsif($orient eq 'W'){$vx1=$x2;$vy1=$y1;$vx2=$x2;$vy2=$y2;$vx3=$x2-$p;$vy3=$ymid;}
  elsif($orient eq 'E'){$vx1=$x1;$vy1=$y1;$vx2=$x1;$vy2=$y2;$vx3=$x1+$p;$vy3=$ymid;}

  # now draw the triangle
  $gd->line($vx1,$vy1,$vx2,$vy2,$fg);
  $gd->line($vx2,$vy2,$vx3,$vy3,$fg);
  $gd->line($vx3,$vy3,$vx1,$vy1,$fg);

  if (my $c = $self->bgcolor){
    $gd->fillToBorder($xmid,$ymid,$fg,$c) if $orient eq 'S' || $orient eq 'N';
    $gd->fillToBorder($x1+1,$ymid,$fg,$c) if $orient eq 'E';
    $gd->fillToBorder($x2-1,$ymid,$fg,$c) if $orient eq 'W';
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::triangle - The "triangle" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws an equilateral triangle when -point is defined.
It draws an isoceles triangle otherwise.  It is possible to draw
the triangle with the base on the N, S, E, or W side.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -point      If true, the triangle         0
              will drawn at the center
              of the range, and not scaled
              to the feature width.

  -orient     On which side shall the       S
              base be? (NSEW)

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
