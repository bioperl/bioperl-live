package Bio::Graphics::Glyph::track;

use strict;
use base qw(Bio::Graphics::Glyph);

# track sets connector to empty
sub connector {
  my $self = shift;
  return $self->SUPER::connector(@_) if $self->all_callbacks;
  return 'none';
}

sub draw {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;

  # the clipping code here prevents poorly-behaving glyphs from
  # drawing outside the track
  my @clip;
  if ($gd->can('clip')) {
    @clip = $gd->clip();
    # glyphs are allowed a slop area of ~3 on either side and 6 on the top and bottom
    # in order to spill out over their boundaries.  Beyond this they start overlapping
    # with other glyphs in an ugly way.
    my @cliprect = ($left-$self->panel->pad_left,
		    $top-6,
		    $self->panel->right+$self->panel->pad_right,
		    $top+$self->layout_height+6);
    $gd->clip(@cliprect);
  }

  my @parts = $self->parts;
  for (my $i=0; $i<@parts; $i++) {
    $parts[$i]->draw_highlight($gd,$left,$top);
    $parts[$i]->draw($gd,$left,$top,0,1);
  }

  $gd->clip(@clip) if @clip;
}

# do nothing for components
# sub draw_component { }

1;


__END__

=head1 NAME

Bio::Graphics::Glyph::track - The "track" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used internally by Bio::Graphics::Panel for laying out
tracks.  It should not be used explicitly.

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

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
