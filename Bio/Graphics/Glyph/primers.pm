package Bio::Graphics::Glyph::primers;
# package to use for drawing something that looks like
# primer pairs.

use strict;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';
use Bio::Graphics::Glyph::generic;

use constant HEIGHT => 4;

# we do not need the default amount of room
#sub calculate_height {
#  my $self = shift;
#  return $self->option('label') ? HEIGHT + $self->labelheight + 2 : HEIGHT;
#}

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
  my $trunc_left  = $x1 < $self->panel->left;
  my $trunc_right = $x2 > $self->panel->right;

  unless ($trunc_left) { 
    $gd->line($x1,$center,$x1 + HEIGHT,$center,$fg);
    $gd->line($x1 + HEIGHT,$center,$x1 + HEIGHT - $a2,$center-$a2,$fg);
    $gd->line($x1 + HEIGHT,$center,$x1 + HEIGHT - $a2,$center+$a2,$fg);
  }

  unless ($trunc_right) {
    $gd->line($x2,$center,$x2 - HEIGHT,$center,$fg);
    $gd->line($x2 - HEIGHT,$center,$x2 - HEIGHT + $a2,$center+$a2,$fg);
    $gd->line($x2 - HEIGHT,$center,$x2 - HEIGHT + $a2,$center-$a2,$fg);
  }

  # connect the dots if requested
  if ($self->connect) {
    my $c = $self->color('connect_color') || $self->bgcolor;
    $gd->line($x1 + ($trunc_left  ? 0 : HEIGHT + 2),$center,
	      $x2 - ($trunc_right ? 0 : HEIGHT + 2),$center,
	      $c);
  }

  # add a label if requested
  $self->draw_label($gd,@_)       if $self->option('label');
  $self->draw_description($gd,@_) if $self->option('description');

}

sub connect {
  my $self = shift;
  return $self->option('connect') if defined $self->option('connect');
  1;  # default
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

  -connect    Whether to connect the      true
              two arrowheads by a line.

  -connect_color  The color to use for the    bgcolor
              connecting line.

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

Allen Day E<lt>day@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
