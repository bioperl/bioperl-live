package Bio::Graphics::Glyph::toomany;
# DAS-compatible package to use for drawing a box

use strict;
use base qw(Bio::Graphics::Glyph::generic);

# draw the thing onto a canvas
# this definitely gets overridden
sub draw {
  my $self = shift;
  my $gd   = shift;
  my ($left,$top) = @_;

  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries($left,$top);

#  $self->filled_oval($gd,$x1,$y1,$x2,$y2);
  for (my $m = 3;$m > 0;$m--){
    my $stack = $m * $self->height / 2;
    $self->unfilled_box($gd,$x1-$stack,$y1-$stack,$x2-$stack,$y2-$stack);
  }

  # add a label if requested
  $self->draw_label($gd,$left,($top-($self->height*1.1))) if $self->option('label');
}

sub label {
  return "too many to display";
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::toomany - The "too many to show" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is intended for features that are too dense to show
properly.  Mostly a placeholder, it currently shows a filled oval.  If
you choose a bump of 0, the ovals will overlap, to give a cloud
effect.

=head2 OPTIONS

There are no glyph-specific options.

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
