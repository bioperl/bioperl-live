package Bio::Graphics::Glyph::dot;
# DAS-compatible package to use for drawing a ring or filled circle

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::generic;
@ISA = 'Bio::Graphics::Glyph::generic';
use constant PI => 3.14159;

sub draw {
  my $self = shift;
#  $self->SUPER::draw(@_);
  my $gd = shift;
  my $fg = $self->fgcolor;

  # now draw a circle
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  my $xmid   = (($x1+$x2)/2);  my $width  = abs($x2-$x1);
  my $ymid   = (($y1+$y2)/2);  my $height = abs($y2-$y1);

  #only point ovals allowed now
  my $r = $self->height ;
    $gd->arc($xmid,$ymid,$r,$r,0,360,$fg);


  if ($self->option('bgcolor')){
    my $c = $self->color('bgcolor');
    $gd->fill($xmid,$ymid,$c);
  }

  #how about a fuse for the bomb?
  #work in degrees, not radians.  So we define PI above
  if(defined $self->option('stem')){
    my $angle = $self->option('stem');

    $gd->line($xmid+($r/PI*sin($angle*PI/180)),
	      $ymid+($r/PI*cos($angle*PI/180)),
	      $xmid+($r*sin($angle*PI/180)),
	      $ymid+($r*cos($angle*PI/180)),$fg);
  }

  $self->draw_label($gd,@_) if $self->option('label');
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::dot - The "dot" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws an ellipse the width of the scaled feature passed,
and height a possibly configured height (See Bio::Graphics::Glyph).

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

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description                  Default
  ------      -----------                  -------

  -point      Whether to draw an ellipse   feature width
              the scaled width of the
              feature or with radius
              point.

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
