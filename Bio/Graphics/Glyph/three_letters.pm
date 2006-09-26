package Bio::Graphics::Glyph::three_letters;
# DAS-compatible package to use for drawing a line of groups of three letters

# $Id$
# Non object-oriented utilities used here-and-there in Bio::Graphics modules

=head1 NAME

Bio::Graphics::Glyph::three_letters - DAS-compatible package to use for drawing a line of groups of three letters

=cut

use strict;
use base qw(Bio::Graphics::Glyph::repeating_shape);

sub pad_top {
  my $self = shift;
  my $top = $self->SUPER::pad_top;
  my $extra = 0.2 * $self->font->height;
  return $top + $extra;
}

sub default_interval
{
  return 20;  
}

sub default_text
{
	return "CAG";
}

sub draw_repeating_shape
{
  my ($self, $gd, $x1, $y1, $x2, $y2, $fg) = @_;
  
  my $text = defined $self->option('text') ? $self->option('text') : $self->default_text();
  
  while (length $text < 3)
  {
    $text .= " ";  
  }
  
  $text = substr($text,0,3);
  my @letters = split //, $text;  
  
  my $oneThird = ($x2-$x1) / 3;
  my $secondLetterX = $x1 + $oneThird;
  my $thirdLetterX = $x1 + 2*$oneThird;

  my $font = $self->option('labelfont') || $self->font;
  $gd->string($font, $x1, $y2-$font->height, $letters[0], $self->fontcolor);
  $gd->string($font, $secondLetterX, $y2-1.7*$font->height, $letters[1], $self->fontcolor);
  $gd->string($font, $thirdLetterX, $y2-$font->height, $letters[2], $self->fontcolor);  
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::three-letters - The "three letters" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws groups of three letters separated by horizontal lines.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------
  -text       The three letters to show     "CAG"

  -width      Width of one letter group     20

  -interval   Interval between              10
              letter groups

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
