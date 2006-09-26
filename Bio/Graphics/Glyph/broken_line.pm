package Bio::Graphics::Glyph::broken_line;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_draw_beak
{
  return 1;  
}

sub default_shear
{
  return 5;  
}

sub default_shear_up
{
  return 1;  
}

sub default_break
{
  return 8;  
}

sub default_extend
{
  return 1;  
}

sub default_size
{
  return 30;  
}

sub default_omit_left
{
  return 0;  
}

sub default_omit_right
{
  return 0;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;
  
  my $shear = defined $self->option('shear') ? $self->option('shear') : $self->default_shear();
  my $shear_up = defined $self->option('shear_up') ? $self->option('shear_up') : $self->default_shear_up();
  my $break = defined $self->option('break') ? $self->option('break') : $self->default_break();
  my $draw_beak = defined $self->option('draw_beak') ? $self->option('draw_beak') : $self->default_draw_beak();
  my $extend= defined $self->option('extend') ? $self->option('extend') : $self->default_extend();
  my $size = defined $self->option('size') ? $self->option('size') : $self->default_size();
  my $omit_left = defined $self->option('omit_left') ? $self->option('omit_left') : $self->default_omit_left();
  my $omit_right = defined $self->option('omit_right') ? $self->option('omit_right') : $self->default_omit_right();

  my $midY = ($y1+$y2)/2;
  
  if ($x2-$x1 < $size)
  {
    $gd->line($x1, $midY, $x2, $midY, $bg);
    return;
  }
  
  my $midX = ($x1+$x2)/2;
  
  my $break_start = $midX - $break/2;
  my $break_end = $midX + $break/2;
  
  my ($x11, $x12, $x21, $x22);
  $x12 = $break_start;
  $x21 = $break_end;

  if ($omit_left)
  {
    $break_start = $x1;
    $break_end = $x1+$break;
    $x21 = $break_end;
    $x22 = ($extend ? $x2 : $x21 + $size - $break);
  }
  elsif ($omit_right)
  {
    $x11 = $x1;
    $x12 = ($extend ? $x2 - $break : $x11 + $size - $break);
    $break_end = $x12+$break;
    $break_start = $x12;
  }
  else
  {
    if ($extend)
    {
      $x11 = $x1;
      $x22 = $x2;
    }
    else
    {
      $x11 = $break_start - ($size - $break) / 2;
      $x22 = $break_end + ($size - $break) / 2;
    }
  }
  
  unless ($omit_left)
  {
    $gd->line($x11, $midY, $x12, $midY, $bg);
  }
  
  my $shear_y = ($shear_up ? $midY - $shear : $midY + $shear);
  $gd->line($break_start, $shear_y, $break_end, $shear_y, $fg);
  if ($draw_beak)
  {
    $midX = ($break_start + $break_end) / 2;
    
    my $beak_y1 = $shear_up ? $midY + $shear/2 : $midY - $shear/2;
    my $beak_y2 = $shear_up ? $midY - $shear/2 : $midY + $shear/2;
    
    $gd->line($midX, $beak_y1, $midX-$shear, $beak_y2, $fg);  
    $gd->line($midX, $beak_y1, $midX+$shear, $beak_y2, $fg);  
  }
  
  unless ($omit_right)
  {
    $gd->line($x21, $midY, $x22, $midY, $bg);
  }
  
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::broken_line - The "broken line" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a straight line whose segment is shifted ('sheared')
up or down. There can be an optional "beak' (two diagonal lines
passing between the main line and its segment).
Either the left or the right side of the main line can be absent.
The line can be of fixed size or extend to take up all available space.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -draw_beak Whether to draw the 'beak'.        1

  -shear    Vertical distance between     	5
			the main line and the segment

  -shear_up Whether to shift the segment 	1
				up or down (1 or 0)

  -break 	Width of the break in the line 	8

  -extend  	Whether to extend the line or   1 
			to keep the length fixed (1 or 0) 

  -size  	Total length of the line and   30 
			the break, if extend is 0

  -omit_left	Whether to omit the left	0
			half of the main line

  -omit_right	Whether to omit the right	0
			half of the main line

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
