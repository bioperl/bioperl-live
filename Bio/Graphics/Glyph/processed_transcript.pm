package Bio::Graphics::Glyph::processed_transcript;

# $Id$

use strict;
use Bio::Graphics::Glyph::transcript2;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::transcript2';
use constant DEFAULT_UTR_COLOR => '#D0D0D0';

# This hack preprocesses the glyph to remove overlaps between UTRs and
# exons.  The exons are clipped so that UTRs have precedence
sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);

  if ($self->option('adjust_exons')) {
    # find the utrs
    my @utrs  = grep {$_->feature->type =~ /utr/i} $self->parts;
    my %seen  = map {$_=>1} @utrs;
    my @other = grep {!$seen{$_}} $self->parts;
    my @clipped_parts;

    for my $p (@other) {
      my $p_right = $p->{left}+$p->{width};
      for my $u (@utrs) {
	my $u_right = $u->{left}+$u->{width};
	if ($u_right > $p->{left} && $u->{left} < $p_right) { # overlaps - clip
	  if ($p->{left} >= $u->{left}) {
	    $p->{left}    = $u_right;
	    $p->{width}   = $p_right - $p->{left};
	  }
	  if ($p_right <= $u_right) {
	    my $new_right = $u->{left};
	    $p->{width}   = $u->{left} - $p->{left};
	  }
	}
      }
    }
    $self->{parts} = [grep {$_->{width}>1} sort {$a->{left}<=>$b->{left}} $self->parts];
  }

  $self;
}

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

sub connector {
  my $self = shift;
  return 'quill' if $self->option('decorate_introns');
  return $self->SUPER::connector(@_);
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

  -connector    Connector type                 undef (false)

  -connector_color
                Connector color                black

  -label        Whether to draw a label	       undef (false)

  -description  Whether to draw a description  undef (false)

  -strand_arrow Whether to indicate            undef (false)
                 strandedness

  -hilite       Highlight color                undef (no color)

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option         Description                  Default
  ------         -----------                  -------

  -thin_utr      Flag.  If true, UTRs will      undef (false)
                 be drawn at 2/3 of the
                 height of CDS segments.

  -utr_color     Color of UTR segments.         Gray #D0D0D0

  -decorate_introns
                 Draw strand with little arrows undef (false)
                 on the intron.

  -adjust_exons  Fix exons so that they don't   undef (false)
                 overlap UTRs

The B<-adjust_exons> option is needed to handle features in which the
exons (SO type "exon") overlaps with the UTRs (SO types
"five_prime_utr" and "three_prime_utr").  The exon parts of the glyph
will be clipped so that it doesn't overlap with the UTR parts.

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
