package Bio::Graphics::Glyph::primers;
# package to use for drawing something that looks like
# primer pairs.

use strict;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';

use constant HEIGHT => 4;

# we do not need the default amount of room
sub calculate_height {
  my $self = shift;
  return $self->option('label') ? HEIGHT + $self->labelheight + 2 : HEIGHT;
}

# override draw method
sub draw {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $fg = $self->fgcolor;
  my $a2 = HEIGHT/2;
  my $center = $y1 + $a2;

  # just draw us as a solid line -- very simple
  if ($x2-$x1 < HEIGHT*2) {
    $gd->line($x1,$center,$x2,$center,$fg);
    return;
  }

  # otherwise draw two pairs of arrows
  # -->   <--
  $gd->line($x1,$center,$x1 + HEIGHT,$center,$fg);
  $gd->line($x1 + HEIGHT,$center,$x1 + HEIGHT - $a2,$center-$a2,$fg);
  $gd->line($x1 + HEIGHT,$center,$x1 + HEIGHT - $a2,$center+$a2,$fg);

  $gd->line($x2,$center,$x2 - HEIGHT,$center,$fg);
  $gd->line($x2 - HEIGHT,$center,$x2 - HEIGHT + $a2,$center+$a2,$fg);
  $gd->line($x2 - HEIGHT,$center,$x2 - HEIGHT + $a2,$center-$a2,$fg);

  # connect the dots if requested
  if ($self->option('connect')) {
    my $c = $self->color('connect_color') || $self->bgcolor;
    $gd->line($x1 + HEIGHT + 2,$center,$x2 - HEIGHT - 2,$center,$c);
  }

  # add a label if requested
  $self->draw_label($gd,@_) if $self->option('label');

}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::primers - The "STS primers" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws two arrows oriented towards each other and connected
by a line of a contrasting color.  The length of the arrows is
immaterial, but the length of the glyph itself corresponds to the
length of the scaled feature.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -connect    Whether to connect the      false
              two arrowheads by a line.

  -connect_color  The color to use for the    bgcolor
              connecting line.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Ace::Sequence>, L<Ace::Sequence::Feature>, L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>, L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,

=head1 AUTHOR

Allen Day <day@cshl.org>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
