package Bio::Graphics::Glyph::transcript2;

use strict;
use Bio::Graphics::Glyph::transcript;
use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::transcript';

use constant MIN_WIDTH_FOR_ARROW => 8;

sub pad_left  {
  my $self = shift;
  my $pad = $self->Bio::Graphics::Glyph::generic::pad_left;
  return $pad unless $self->feature->strand < 0;
  my $first = ($self->parts)[0] or return $pad;
  my @rect  = $first->bounds();
  my $width = abs($rect[2] - $rect[0]);
  return $self->SUPER::pad_left if $width < MIN_WIDTH_FOR_ARROW;
  return 0;
}

sub pad_right  {
  my $self = shift;
  my $pad = $self->Bio::Graphics::Glyph::generic::pad_right;
  my $last = ($self->parts)[-1] or return $pad;
  my @rect  = $last->bounds();
  my $width = abs($rect[2] - $rect[0]);
  return $self->SUPER::pad_right if $width < MIN_WIDTH_FOR_ARROW;
  return $pad;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my ($left,$top) = @_;
  my @rect = $self->bounds(@_);

  my $width = abs($rect[2] - $rect[0]);
  my $filled = defined($self->{partno}) && $width >= MIN_WIDTH_FOR_ARROW;

  if ($filled) {
    my $f = $self->feature;

    if ($f->strand < 0 
	&& (!$self->is_recursive
	 || $self->{partno} == 0)) { # first exon, minus strand transcript
      $self->filled_arrow($gd,-1,@rect);
      $self->{filled}++;
    } elsif ($f->strand >= 0 
	     && (!$self->is_recursive
		 || $self->{partno} == $self->{total_parts}-1)) { # last exon, plus strand
      $self->filled_arrow($gd,+1,@rect);
      $self->{filled}++;
    } else {
      $self->SUPER::draw_component($gd,@_);
    }
  }

  else {
    $self->SUPER::draw_component($gd,@_);
  }

}

# override option() for force the "hat" type of connector
sub connector {
  return 'hat';
}

sub draw_connectors {
  my $self = shift;
  my @parts = $self->parts;
  if ($self->{filled} || $parts[0]->{filled} || $parts[-1]->{filled}) {
    $self->Bio::Graphics::Glyph::generic::draw_connectors(@_);
  } else {
    $self->SUPER::draw_connectors(@_);
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
L<Bio::Graphics::Track>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,

=head1 AUTHOR

Lincoln Stein <lstein@cshl.org>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
