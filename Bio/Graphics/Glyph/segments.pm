package Bio::Graphics::Glyph::segments;

use strict;
use Bio::Location::Simple;
use Bio::Graphics::Glyph::generic;
use Bio::Graphics::Glyph::segmented_keyglyph;
use vars '@ISA','$VERSION';
@ISA = qw( Bio::Graphics::Glyph::segmented_keyglyph
	   Bio::Graphics::Glyph::generic
	 );
$VERSION = '1.00';

# group sets connector to 'solid'
sub connector {
  my $self = shift;
  return $self->SUPER::connector(@_) if $self->all_callbacks;
  return $self->SUPER::connector(@_) || 'solid';
}
# never allow our components to bump
sub bump {
  my $self = shift;
  return $self->SUPER::bump(@_) if $self->all_callbacks;
  return 0;
}
sub label {
  my $self = shift;
  return $self->SUPER::label(@_) if $self->all_callbacks;
  return unless $self->{level} == 0;
  return $self->SUPER::label(@_);
}
sub description {
  my $self = shift;
  return $self->SUPER::description(@_) if $self->all_callbacks;
  return unless $self->{level} == 0;
  return $self->SUPER::description(@_);
}

# Override _subseq() method to make it appear that a top-level feature that
# has no subfeatures appears as a feature that has a single subfeature.
# Otherwise at high mags gaps will be drawn as components rather than
# as connectors.  Because of differing representations of split features
# in Bio::DB::GFF::Feature and Bio::SeqFeature::Generic, there is
# some breakage of encapsulation here.
sub _subseq {
  my $self    = shift;
  my $feature = shift;
  my @subseq  = $self->SUPER::_subseq($feature);
  return @subseq if @subseq;
  if ($self->level == 0 && !@subseq && !eval{$feature->compound}) {
    my($start,$end) = ($feature->start,$feature->end);
    ($start,$end) = ($end,$start) if $start > $end; # to keep Bio::Location::Simple from bitching
    return Bio::Location::Simple->new(-start=>$start,-end=>$end);
  } else {
    return;
  }
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::segments - The "segments" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing features that consist of discontinuous
segments.  Unlike "graded_segments" or "alignment", the segments are a
uniform color and not dependent on the score of the segment.

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

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

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
