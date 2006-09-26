package Bio::Graphics::Glyph::dashed_line;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_linewidth
{
  return 1;  
}

sub default_dash_size
{
  return 6;  
}

sub default_space_size
{
  return 3;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $fg = $self->fgcolor;
  
  my $midY = ($y1+$y2) / 2;
  
  my $linewidth = defined $self->option('linewidth') ? $self->option('linewidth')  : $self->default_linewidth();
  my $dash_size = defined $self->option('dash_size') ? $self->option('dash_size') : $self->default_dash_size();
  my $space_size = defined $self->option('space_size') ? $self->option('space_size') : $self->default_space_size();
  my $space_color = $self->option('space_color');
  my $shear = $self->option('shear') || "";
  $space_color = $self->factory->translate_color($space_color) if $space_color;
  
  my ($x, $_y1, $_y2);
  $x = $x1;
  while ($x<$x2)
  {
    my $newX = $x+$dash_size;
    $newX = $x2 if $newX > $x2;
    if ($shear == 1)
    {
      $_y1 = $midY-$linewidth;
      $_y2 = $midY;
    }
    else
    {
      $_y1 = $midY - $linewidth/2;  
      $_y2 = $midY + $linewidth/2;  
    }
    $self->filled_box($gd,$x,$_y1,$newX,$_y2,$fg,$fg);
    last if $newX >= $x2;
    
    $x = $newX;
    $newX = $x+$space_size;
    $newX = $x2 if $newX > $x2;
    if ($space_color)
    {
      if ($shear == 1)
      {
        $_y1 = $midY;
        $_y2 = $midY+$linewidth;
      }
      else
      {
        $_y1 = $midY - $linewidth/2;  
        $_y2 = $midY + $linewidth/2;  
      }
      $self->filled_box($gd, $x,$_y1,$newX,$_y2,$space_color,$space_color);
    }
    $x = $newX;
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::dashed_line - The "dashed line" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a dashed line. The lengths of the dash and the space are configurable.
The space can be filled with a different color, thus making a two-colored line.
Also, the two colors can be "sheared".

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -dash_size  Width of one dash              6

  -space_size Width of one space             3
              between dashes       

  -space_color Color of the space            none 
              between dashes       

  -shear      Whether to use shearing       0
              (1 or 0)

  -linewidth  Standard option, but          1
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
