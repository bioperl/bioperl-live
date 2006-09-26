package Bio::Graphics::Glyph::pinsertion;
# package to use for drawing P insertion as a triangle
# p insertion is a point (one base).

use strict;
use GD;
use base qw(Bio::Graphics::Glyph::generic);

sub box {
  my $self = shift;
  my $half = $self->insertion_width/2;
  return ($self->left-$half,$self->top,$self->right+$half,$self->bottom);
}

sub insertion_width {
  my $self = shift;
  return $self->option('insertion_width') || 6;
}

# override draw method
sub draw {
  my $self = shift;

  my $gd = shift;
  my ($left,$top) = @_;

  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $height = $self->height;

  my $half = $self->insertion_width/2;
  my $fill = $self->bgcolor;
  my $border = $self->fgcolor;

  my $poly_pkg = $self->polygon_package;
  my $poly     = $poly_pkg->new();
  if ($self->feature->strand > 0) { #plus strand
      $poly->addPt($x1 - $half, $y1);
      $poly->addPt($x1 + ($half), $y1);
      $poly->addPt($x1, $y2); #pointer down
  } else {
      $poly->addPt($x1, $y1); #pointer up
      $poly->addPt($x1 - $half, $y2);
      $poly->addPt($x1 + ($half), $y2);
  }
  $gd->filledPolygon($poly,$fill);
  $gd->polygon($poly,$border);

  # add a label if requested
  $self->draw_label($gd,$left,$top)       if $self->option('label');
  $self->draw_description($gd,$left,$top) if $self->option('description');
}


1;

=head1 NAME

Bio::Graphics::Glyph::pinsertion - The "Drosophila P-element Insertion" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph was designed to show P-element insertions in the Drosophila
genome, but in fact is suitable for any type of zero-width feature.
Also see the triangle glyph.

=head2 OPTIONS

In addition to the generic options, this glyph recognizes:

 Option Name       Description              Default
 -----------       -----------              -------

 -insertion_width  Width of glyph in pixels    3

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

Allen Day E<lt>day@cshl.orgE<gt>, Shengqiang Shu E<lt>sshu@bdgp.lbl.govE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory, BDGP

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
