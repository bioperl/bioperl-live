package Bio::Graphics::Glyph::transcript2;

# $Id$

use strict;
use base qw(Bio::Graphics::Glyph::transcript);

use constant MIN_WIDTH_FOR_ARROW => 8;

sub extra_arrow_length {
  my $self = shift;
#  return 0 unless $self->{level} == 0;
  my $strand = $self->feature->strand || 0;
  $strand *= -1 if $self->{flip};
  return 0 unless $strand < 0;
  my $first = ($self->parts)[0] || $self;
  my @rect  = $first->bounds();
  my $width = abs($rect[2] - $rect[0]);
  return 0 if $width >= MIN_WIDTH_FOR_ARROW;
  return $self->arrow_length;
}

sub pad_left  {
   my $self = shift;
#   return 0 unless $self->{level} == 0;
   my $pad = $self->Bio::Graphics::Glyph::generic::pad_left;
   my $extra_arrow_length = $self->extra_arrow_length;
   if ($self->label_position eq 'left' && $self->label) {
     return $extra_arrow_length+$pad;
   } else {
     return $extra_arrow_length > $pad ? $extra_arrow_length : $pad;
   }
}

sub pad_right  {
  my $self = shift;
#  return 0 unless $self->{level} == 0;
  my $pad = $self->Bio::Graphics::Glyph::generic::pad_right;
  return $pad if $self->{level} > 0;
  my $last = ($self->parts)[-1] || $self;
  my @rect  = $last->bounds();
  my $width = abs($rect[2] - $rect[0]);
  return $self->SUPER::pad_right if $width < MIN_WIDTH_FOR_ARROW;
  return $pad
}

sub draw_connectors {
  my $self = shift;
  my ($gd,$dx,$dy) = @_;

  my $part;
  my $strand = $self->feature->strand;
  $strand   *= -1 if $self->{flip};  #sigh
  if (my @parts  = $self->parts) {
    $part   = $strand >= 0 ? $parts[-1] : $parts[0];
  } elsif ($self->feature_has_subparts) {
    # no parts -- so draw an intron spanning whole thing
    my($x1,$y1,$x2,$y2) = $self->bounds(0,0);
    $self->_connector($gd,$dx,$dy,$x1,$y1,$x1,$y2,$x2,$y1,$x2,$y2);
    $part = $self;
  } else {
    return;
  }
  my @rect   = $part->bounds();
  my $width  = abs($rect[2] - $rect[0]);
  my $filled = $width >= MIN_WIDTH_FOR_ARROW;

  if ($filled) {
    $self->Bio::Graphics::Glyph::generic::draw_connectors(@_);
  } else {
    $self->SUPER::draw_connectors(@_);
  }
}

sub draw_component {
  my $self = shift;
  return unless $self->level > 0;

  my $gd = shift;
  my ($left,$top) = @_;
  my @rect = $self->bounds(@_);

  my $f      = $self->feature;
  my $strand = $f->strand;
  my $str    = $strand * ($self->{flip} ? -1 : 1);

  my $width = abs($rect[2] - $rect[0]);
  my $filled = defined($self->{partno}) && $width >= MIN_WIDTH_FOR_ARROW;
  my ($pwidth) = $gd->getBounds;
  $filled = 0 if $str < 0 && $rect[0] < $self->panel->pad_left;
  $filled = 0 if $str > 0 && $rect[2] > $pwidth - $self->panel->pad_right;

  if ($filled) {
    my ($first,$last)  = ($self->{partno} == 0 , $self->{partno} == $self->{total_parts}-1);
    ($first,$last)     = ($last,$first) if $self->{flip};

    if ($strand < 0 && $first) { # first exon, minus strand transcript
      $self->filled_arrow($gd,-1,@rect);
    } elsif ($strand >= 0 && $last) { # last exon, plus strand
      $self->filled_arrow($gd,+1,@rect);
    } else {
      $self->SUPER::draw_component($gd,@_);
    }
  }

  else {
    $self->SUPER::draw_component($gd,@_);
  }

}

sub bump {
  my $self = shift;
  return $self->SUPER::bump(@_) if $self->all_callbacks;
  return 0;  # never allow our components to bump
}

1;


__END__

=head1 NAME

Bio::Graphics::Glyph::transcript2 - The "transcript2" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing transcripts.  It is like "transcript"
except that if there is sufficient room the terminal exon is shaped
like an arrow in order to indicate the direction of transcription.  If
there isn't enough room, a small arrow is drawn.

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
