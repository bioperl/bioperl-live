package Bio::Graphics::Glyph::repeating_shape;
# DAS-compatible package to use for drawing a line of repeating shapes

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_width
{
  return 10;  
}

sub default_interval
{
  return 10;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  
  my $width = defined $self->option('width') ? $self->option('width') : $self->default_width;
  my $interval = defined $self->option('interval') ? $self->option('interval') : $self->default_interval;

  # find the center and vertices
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $bWidth = $x2-$x1;
  
  if ($bWidth < $width)
  {
    $self->draw_repeating_shape($gd,$x1,$y1,$x2,$y2,$fg);
    return;
  }
  
  if ($bWidth < $width+2*$interval)
  {
    my $leftoverInterval = $bWidth - $width;
    my $halfInt = $leftoverInterval/2;
    $halfInt = 0 unless $interval;
    
    $gd->line($x1,$y2,$x1+$halfInt,$y2,$fg);
    $self->draw_repeating_shape($gd,$x1+$halfInt,$y1,$x2-$halfInt,$y2,$fg);
    $gd->line($x2-$halfInt,$y2,$x2,$y2,$fg);
    return;
  }
  
  my $count = int ($bWidth / ($width+$interval));
  my $leftoverInterval = $bWidth % ($width+$interval)+$interval;
  
  my $halfInt = $leftoverInterval/2;
  $halfInt = 0 unless $interval;
  $gd->line($x1,$y2,$x1+$halfInt,$y2,$fg);
  foreach (my $i=1; $i<=$count; $i++)
  {
    my $shapeStart = $x1 + $halfInt + ($i-1)*($width+$interval);
    $self->draw_repeating_shape($gd,$shapeStart,$y1,$shapeStart+$width,$y2,$fg);
    if ($i < $count)
    {
      $gd->line($shapeStart+$width,$y2,$shapeStart+$width+$interval,$y2,$fg);  
    }
  }
  $gd->line($x2-$halfInt,$y2,$x2,$y2,$fg);
}

sub draw_repeating_shape
{
	warn "Subclasses must implement 'draw_repeating_shape'!\n";
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::repeating_shape - A glyph that draws the same shape repeatedly. 

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is a generic superclass for drawing the same shape repeatedly.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -width      Width of one tooth            10

  -interval   Interval between teeth        10

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
