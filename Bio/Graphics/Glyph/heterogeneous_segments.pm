package Bio::Graphics::Glyph::heterogeneous_segments;

# this glyph acts like graded_segments but the bgcolor of each segment is
# controlled by the source field of the feature. Use the source field name
# to set the background color:
# -waba_strong => 'blue'
# -waba_weak   => 'red'
# -waba_coding => 'green' 

use strict;
use Bio::Graphics::Glyph::graded_segments;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::graded_segments';

# override draw method to calculate the min and max values for the components
sub draw {
  my $self = shift;

  # bail out if this isn't the right kind of feature
  # handle both das-style and Bio::SeqFeatureI style,
  # which use different names for subparts.
  my @parts = $self->parts;
  return $self->SUPER::draw(@_) unless @parts;

  # figure out the colors
  $self->{source2color} ||= {};
  my $fill = $self->bgcolor;
  for my $part (@parts) {
    my $s = eval { $part->feature->source_tag } or next;
    $self->{source2color}{$s} ||= $self->color(lc($s)."_color") || $fill;
    $part->{partcolor} = $self->{source2color}{$s};
  }

  $self->Bio::Graphics::Glyph::generic::draw(@_);
}


# synthesize a key glyph
sub keyglyph {
  my $self = shift;
  
  my $scale = 1/$self->scale;  # base pairs/pixel

  # two segments, at pixels 0->50, 60->80
  my $offset = $self->panel->offset;

  my $feature =
    Bio::Graphics::Feature->new(
				-segments=>[ [ 0*$scale +$offset,25*$scale+$offset],
					     [ 25*$scale +$offset,50*$scale+$offset],
					     [ 50*$scale+$offset, 75*$scale+$offset]
					   ],
				-name => $self->option('key'),
				-strand => '+1');
  my @sources = grep {/_color$/} $self->factory->options;
  foreach (@sources) {s/_color$//}
  ($feature->segments)[0]->source_tag($sources[1]);
  ($feature->segments)[1]->source_tag($sources[0]);
  ($feature->segments)[2]->source_tag($sources[2]);
  my $factory = $self->factory->clone;
  $factory->set_option(label => 1);
  $factory->set_option(bump  => 0);
  $factory->set_option(connector  => 'solid');
  my $glyph = $factory->make_glyph(0,$feature);
  return $glyph;
}

# component draws a shaded box
sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($left,$top) = @_;
  my $color = $self->{partcolor};
  my @rect = $self->bounds(@_);
  $self->filled_box($gd,@rect,$color,$color);
}

1;

=head1 NAME

Bio::Graphics::Glyph::heterogeneous_segments - The "heterogeneous_segments" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph acts like graded_segments but the bgcolor of each segment (sub-feature)
can be individually set using the source field of the feature.

Each segment type color is specified using the following nomenclature:

 -{source}_color => $color

For example, if the feature consists of a gene containing both
confirmed and unconfirmed exons, you can make the confirmed exons
green and the unconfirmed ones red this way:

  -confirmed_color   => 'green',
  -unconfirmed_color => 'red'

=head2 OPTIONS

The following options are standard among all Glyphs.  See
L<Bio::Graphics::Glyph> for a full explanation.

  Option      Description                      Default
  ------      -----------                      -------

  -fgcolor      Foreground color	       black

  -outlinecolor	Synonym for -fgcolor

  -bgcolor      Background color               turquoise

  -fillcolor    Synonym for -bgcolor

  -linewidth    Line width                     1

  -height       Height of glyph		       10

  -font         Glyph font		       gdSmallFont

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

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
