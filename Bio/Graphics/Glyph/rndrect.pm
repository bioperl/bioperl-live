package Bio::Graphics::Glyph::rndrect;

use strict;
use base 'Bio::Graphics::Glyph::generic';

# override draw_component to draw an round edge rect rather than a rectangle
sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($left,$top) = @_;
  my($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);#$self->bounds(@_);

  my $poly_pkg = $self->polygon_package;
  my $poly     = $poly_pkg->new();
  my $boxheight = $y2 - $y1;

  if (($x2-$x1) > 3) {
      $poly->addPt($x1+1, $y1+1);
      $poly->addPt($x1+2, $y1);
      $poly->addPt($x2-2, $y1);
      $poly->addPt($x2-1, $y1+1);
      $poly->addPt($x2, $y1 + $boxheight / 2)
        if (($y2 - $y1) > 6);

      $poly->addPt($x2-1, $y2-1);
      $poly->addPt($x2-2, $y2);
      $poly->addPt($x1+2, $y2);
      $poly->addPt($x1+1, $y2-1);
      $poly->addPt($x1, $y1 + $boxheight / 2)
        if (($y2 - $y1) > 6);
  } else {
      $poly->addPt($x1, $y1);
      $poly->addPt($x2, $y1);

      $poly->addPt($x2, $y2);
      $poly->addPt($x1, $y2);
  }

  $gd->filledPolygon($poly,$self->fillcolor);
  $gd->polygon($poly,$self->fgcolor);
}

# group sets connector to 'solid'
sub connector {
  my $self = shift;
  return $self->SUPER::connector(@_) if $self->all_callbacks;
  return 'solid';
}

sub bump {
  my $self = shift;
  return $self->SUPER::bump(@_) if $self->all_callbacks;
  return 0;
}


1;


=head1 NAME

Bio::Graphics::Glyph::rndrect - The "round rect" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph was designed to show seq features in round edge rectangles.
The glyph will be a rectangle if its width is E<lt> 4 pixels

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

Shengqiang Shu E<lt>sshu@bdgp.lbl.govE<gt>

Copyright (c) 2001 BDGP

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
