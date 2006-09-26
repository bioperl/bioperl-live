package Bio::Graphics::Glyph::transcript;
# $Id$

use strict;
use base qw(Bio::Graphics::Glyph::segments);

sub pad_left  {
  my $self = shift;
  return 0 unless $self->{level} == 0;
  my $pad  = $self->SUPER::pad_left;
  my $strand = $self->feature->strand;
  return $pad unless defined $strand && $strand < 0;
  return $self->arrow_length > $pad ? $self->arrow_length : $pad;
}

sub pad_right {
  my $self = shift;
  return 0 unless $self->{level} == 0;
  my $pad = $self->SUPER::pad_right;
  my $strand = $self->feature->strand;
  return $pad unless defined($strand) && $strand > 0;
  return $self->arrow_length > $pad ? $self->arrow_length : $pad;
}

sub draw_component {
  my $self = shift;
  return unless $self->level > 0;
  $self->SUPER::draw_component(@_);
}

sub part_label_merge {
  my $self = shift;
  my $label = $self->SUPER::part_label_merge;
  return $label if defined $label;
  1;
}

sub draw_connectors {
  my $self = shift;
  my $gd = shift;
  my ($left,$top) = @_;

  $self->SUPER::draw_connectors($gd,$left,$top);
  my @parts = $self->parts; # or return;

  # H'mmm.  No parts.  Must be in an intron, so draw intron
  # spanning entire range
  if (!@parts) {
    return unless $self->feature_has_subparts;
    my($x1,$y1,$x2,$y2) = $self->bounds(0,0);
    $self->_connector($gd,$left,$top,$x1,$y1,$x1,$y2,$x2,$y1,$x2,$y2);
    @parts = ($self);
  }

  # flip argument makes this confusing
  # certainly there's a simpler way to express this idea
  my $strand    = $self->feature->strand;
  my ($first,$last) = ($parts[0],$parts[-1]);

  ($first,$last) = ($last,$first) if exists $self->{flip};

  if ($strand > 0) {
    my($x1,$y1,$x2,$y2) = $last->bounds(@_);
    my $center = ($y2+$y1)/2;
    $self->{flip} ?
	$self->arrow($gd,$x1,$x1-$self->arrow_length,$center)
      :
	$self->arrow($gd,$x2,$x2+$self->arrow_length,$center);
  }

  elsif ($strand < 0) {
    my($x1,$y1,$x2,$y2) = $first->bounds(@_);
    my $center = ($y2+$y1)/2;
    $self->{flip } ?
	$self->arrow($gd,$x2,$x2+$self->arrow_length,$center)
      :
	$self->arrow($gd,$x1,$x1 - $self->arrow_length,$center);
  }
}

sub arrow_length {
  my $self = shift;
  return $self->option('arrow_length') || 8;
}

# override option() for force the "hat" type of connector
sub connector {
  my $self = shift;
  return $self->SUPER::connector(@_) if $self->all_callbacks;
  return ($self->option('connector') || 'hat');
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::transcript - The "transcript" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing transcripts.  It is essentially a
"segments" glyph in which the connecting segments are hats.  The
direction of the transcript is indicated by an arrow attached to the
end of the glyph.

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

  -hilite       Highlight color                undef (no color)

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option         Description                  Default
  ------         -----------                  -------

  -arrow_length  Length of the directional   8
                 arrow.

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
