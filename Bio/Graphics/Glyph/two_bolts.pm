package Bio::Graphics::Glyph::two_bolts;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_bolt_height
{
  return 10;  
}

sub default_bolt_length
{
  return 20;  
}

sub default_remainder_length
{
  return 10;  
}

sub default_bolt_color
{
  return 'red';  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $fg = $self->fgcolor;
  
  my $midY1 = $y1+($y2-$y1) / 3;
  my $midY2 = $y1 + 2*($y2-$y1) / 3;
    
  my $bolt_color = defined $self->option('bolt_color') ? $self->option('bolt_color')  : $self->default_bolt_color();
  $bolt_color = $self->factory->translate_color($bolt_color);

  my $bolt_height = defined $self->option('bolt_height') ? $self->option('bolt_height')  : $self->default_bolt_height();  
  my $bolt_length = defined $self->option('bolt_length') ? $self->option('bolt_length')  : $self->default_bolt_length();  
  my $remainder_length = defined $self->option('remainder_length') ? $self->option('remaindert_length')  : $self->default_remainder_length();  

  if ($x2-$x1 < $bolt_length+$remainder_length)
  {
    $gd->line($x1, $y1, $x2, $y2, $bolt_color);
    return;
  }
  
  my $bolt_start = $x2-$bolt_length-$remainder_length;
  my $step = $bolt_length / 8;
  my $shift = $bolt_height/2;
  $gd->line($x1, $midY1, $bolt_start, $midY1, $fg);
  $self->draw_bolt($gd, $bolt_start, $step, $midY1, $shift, $bolt_color);
  $gd->line($x2-$remainder_length, $midY1, $x2, $midY1, $fg);

  $bolt_start = $x1+$remainder_length;
  $gd->line($x1, $midY2, $bolt_start, $midY2, $fg);
  $self->draw_bolt($gd, $bolt_start, $step, $midY2, $shift, $bolt_color);
  $gd->line($x1+$bolt_length+$remainder_length, $midY2, $x2, $midY2, $fg);
  
  
}

sub draw_bolt
{
  my ($self, $gd, $bolt_start, $step, $y, $shift, $bolt_color) = @_;
  $gd->line($bolt_start, $y, $bolt_start+$step, $y-$shift, $bolt_color);
  $gd->line($bolt_start+$step, $y-$shift, $bolt_start+3*$step, $y+$shift, $bolt_color);
  $gd->line($bolt_start+3*$step, $y+$shift, $bolt_start+5*$step, $y-$shift, $bolt_color);
  $gd->line($bolt_start+5*$step, $y-$shift, $bolt_start+7*$step, $y+$shift, $bolt_color);
  $gd->line($bolt_start+7*$step, $y+$shift, $bolt_start+8*$step, $y, $bolt_color);
  
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::two_bolts - The "two bolts" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws two "bolts" on a line. They look like this;

--------/\/\/\--
--/\/\/\--------

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -bolt_height  Height of the bolt          10

  -bolt_length  Length of the bolt          20

  -bolt_color   Color of the bolt           red

  -remainder_length
                Length of the short line    10

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
