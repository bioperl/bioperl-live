package Bio::Graphics::Glyph::processed_transcript;

# $Id$

use strict;
use Bio::Graphics::Glyph::transcript2;
use vars '@ISA','$VERSION';
@ISA = 'Bio::Graphics::Glyph::transcript2';
$VERSION = '1.0';
use constant DEFAULT_UTR_COLOR => '#D0D0D0';

sub is_utr {
  my $self = shift;
  return $self->feature->primary_tag =~ /UTR|untranslated_region/i;
}

sub thin_utr {
  my $self = shift;
  $self->option('thin_utr');
}

sub utr_color {
  my $self = shift;
  return $self->color('utr_color') if $self->option('utr_color');
  return $self->factory->translate_color(DEFAULT_UTR_COLOR);
}

sub height {
  my $self = shift;
  my $height    = $self->SUPER::height;
  return $height unless $self->thin_utr;
  return $self->is_utr ? int($height/1.5+0.5) : $height;
}

sub pad_top {
  my $self = shift;
  my $pad_top = $self->SUPER::pad_top;
  return $pad_top unless $self->thin_utr && $self->is_utr;
  return $pad_top + int(0.167*$self->SUPER::height + 0.5);
}

sub bgcolor {
  my $self = shift;
  return $self->SUPER::bgcolor unless $self->is_utr;
  return $self->utr_color;
}


1;


__END__

=head1 NAME

Bio::Graphics::Glyph::processed_transcript - The sequence ontology transcript glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing processed transcripts that have both
CDS and UTR segments.  The CDS is drawn in the background color, and
the UTRs are drawn in an alternate color selected by the utr_color
option.  In addition, you can make the UTRs thinner than the CDS by
setting the "thin_utr" option.

For this glyph to produce the desired results, you should pass it a
compound Bio::SeqFeature that has subfeatures of primary_tag "CDS" and
"UTR".  In fact, you may give it more specific types of UTR, including
5'-UTR, 3'-UTR, or the Sequence Ontology terms "untranslated_region,"
"five_prime_untranslated_region," and
"three_prime_untranslated_region."

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

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option         Description                  Default
  ------         -----------                  -------

  -thin_utr      Flag.  If true, UTRs will      0 (false)
                 be drawn at 2/3 of the
                 height of CDS segments.

  -utr_color     Color of UTR segments.         Gray #D0D0D0

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
