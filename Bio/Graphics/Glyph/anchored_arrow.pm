package Bio::Graphics::Glyph::anchored_arrow;
# package to use for drawing an arrow

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::arrow;
@ISA = 'Bio::Graphics::Glyph::arrow';

sub draw_label {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->label or return;
  my $label_align = $self->option('label_align');
  if ($label_align && ($label_align eq 'center' || $label_align eq 'right')) {
      my $x = $self->left + $left;
      my $font = $self->option('labelfont') || $self->font;
      my $middle = $self->left + $left + ($self->right - $self->left) / 2;
      my $label_width = $font->width * length($label);
      if ($label_align eq 'center') {
          my $new_x = $middle - $label_width / 2;
          $x = $new_x if ($new_x > $x);;
      }
      else {
          my $new_x = $left + $self->right - $label_width;
          $x = $new_x if ($new_x > $x);
      }
      $x = $self->panel->left + 1 if $x <= $self->panel->left;
      #detect collision (most likely no bump when want centering label)
      #lay down all features on one line e.g. cyto bands
      return if (!$self->option('bump') && ($label_width + $x) > $self->right);
      $gd->string($font,
                  $x,
                  $self->top + $top,
                  $label,
                  $self->fontcolor);
  }
  else {
      $self->SUPER::draw_label(@_);
  }
}

sub arrowheads {
  my $self = shift;
  my ($ne,$sw,$base_e,$base_w);
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);

  my $gstart  = $x1;
  my $gend    = $x2;
  my $pstart  = $self->panel->left;
  my $pend    = $self->panel->right-1;

  if ($gstart <= $pstart) {  # off left end
    $sw = 1;
  }
  if ($gend >= $pend) { # off right end
    $ne = 1;
  }
  return ($sw,$ne,!$sw,!$ne);
}

sub no_trunc {
  !shift->option('no_arrows');
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::anchored_arrow - The "anchored_arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws an arrowhead which is anchored at one or both ends
(has a vertical base) or has one or more arrowheads.  The arrowheads
indicate that the feature does not end at the edge of the picture, but
continues.  For example:

    |-----------------------------|          both ends in picture
 <----------------------|                    left end off picture
         |---------------------------->      right end off picture
 <------------------------------------>      both ends off picture

You can also set the glyph so that the end is just truncated at the
end of the picture.

         |-----------------------------

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

In addition to the standard options, this glyph recognizes the following:

  Option         Description                Default

  -tick          draw a scale               0 (false)

  -rel_coords    use relative coordinates   0 (false)
                 for scale

  -no_arrows     don't draw an arrow when   0 (false)
                 glyph is partly offscreen

The argument for B<-tick> is an integer between 0 and 2 and has the same
interpretation as the B<-tick> option in Bio::Graphics::Glyph::arrow.

If B<-rel_coords> is set to a true value, then the scale drawn on the
glyph will be in relative (1-based) coordinates relative to the beginning
of the glyph.

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
