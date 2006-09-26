package Bio::Graphics::Glyph::splice_site;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

use constant PWIDTH => 3;

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($left,$top) = @_;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_); 

  my $center = int(0.5+($x1+$x2)/2);
  my $direction = $self->option('direction');
  
  my $height    = $y2 - $y1;
  my $fraction  = $self->option('height_fraction') || 1.0;
  my $bottom    = $y2;
  $top          = $y2 - $fraction * $height;
  
  # draw the base
  my $fgcolor = $self->fgcolor;
  $gd->line($center,$bottom,$center,$top,$fgcolor);

  if ($direction eq 'right') {
    $gd->line($center,$top,$center + PWIDTH,$top,$fgcolor);
  } elsif ($direction eq 'left') {
    $gd->line($center,$top,$center - PWIDTH,$top,$fgcolor);
  }
    
}

1;


=head1 NAME

Bio::Graphics::Glyph::splice_site - The "splice_site" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph was designed to show an inverted "L" representing splice
donors and acceptors.  The vertical part of the L points downwards and
is positioned in the center of the range (even if the range is quite
large).  

In addition to the usual glyph options, this glyph recognizes:

   Option            Value              Description
   ------            -----              -----------

   direction         "left" or "right"  direction the short part of the L
                                        points

   height_fraction   0.0 - 1.0          fractional height of the glyph,
                                        usually a callback that uses the
                                        feature's score to determine its
                                        height

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
L<Bio::Graphics::Glyph::chromosome>,
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

Xiaokang Pan E<lt>pan@cshl.orgE<gt>

Copyright (c) 2001 BDGP

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
