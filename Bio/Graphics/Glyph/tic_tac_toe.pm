package Bio::Graphics::Glyph::tic_tac_toe;

use strict;
use base qw(Bio::Graphics::Glyph::generic);

sub default_mode
{
  return 'x';  
}

sub default_size
{
  return 10;  
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  
  my $fg = $self->fgcolor;
  
  my $size = defined $self->option('size') ? $self->option('size') : $self->default_size();

  my $mode= defined $self->option('mode') ? $self->option('mode') : $self->default_mode();
  
  my $midY = ($y1+$y2)/2;
  
  for (my $i=0; $i<($x2-$x1)/$size; $i++)
  {
    my $start = $x1+$i*$size;
    my $end = $x1+($i+1)*$size;
    if ($mode eq "x" || ($mode eq "xo" && $i%2==0))
    {
      $gd->line($start, $midY-$size/2, $end, $midY+$size/2, $fg);
      $gd->line($end, $midY-$size/2, $start, $midY+$size/2, $fg);
    }
    elsif ($mode eq "o" || ($mode eq "xo" && $i%2==1))
    {
      $gd->ellipse(($start+$end)/2, $midY, $size, $size, $fg);
    }
  }   
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::tic_tac_toe - The "tic-tac-toe" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws a sequence of either 'xxx', 'ooo' or 'xoxo',
depending on the value of 'mode'.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -mode       One of 'x', 'o', or 'xo'.     'x'

  -size       Size of either 'x' or 'o'     10

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
