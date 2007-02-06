package Bio::Graphics::Glyph::fixedwidth;

use strict;
use base 'Bio::Graphics::Glyph::box';
use Carp 'cluck';
use constant TOP_SPACING => 8;

sub pad_left {
  my $self   = shift;
  my $pl     = $self->SUPER::pad_left;
  my $width  = $self->width;
  my $needed = $self->width_needed;
  my $extra  = ($needed-$width)/2 + 1;
  if ($extra > $pl) {
    return $extra;
  } else {
    return $pl;
  }
}

sub pad_right {
  my $self   = shift;
  my $pr     = $self->SUPER::pad_right;
  my $width  = $self->width;
  my $needed = $self->width_needed;
  my $extra  = ($needed-$width)/2;
  $extra = 0 if $extra < 0;
  if ($extra > 0 && $extra > $pr) {
    return $extra;
  } else {
    return $pr-$extra;
  }
}

sub width_needed {
  my $self = shift;
  $self->option('fixed_width');
}

sub height_needed {
  my $self = shift;
  my $h = $self->option('fixed_height');
  return $h if defined $h;
  return $self->SUPER::height;
}

sub height {
  shift->height_needed;
}

sub span_height {
  shift->option('span_height') || 3;
}

sub top_spacing {
  my $self = shift;
  my $spacing = $self->option('fixed_gap');
  return $spacing if defined $spacing;
  return TOP_SPACING;
}

sub pad_top {
  my $self    = shift;
  my $top  = $self->SUPER::pad_top;
  return $top + $self->top_spacing + $self->span_height;
}

sub maxdepth { 0 }

sub draw {
  my $self = shift;
  my $gd       = shift;
  my ($dx,$dy) = @_;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $width        = $self->width_needed;
  my $xmid         = ($x1+$x2) / 2;

  my $top    = $y1 - $self->span_height-$self->top_spacing;
  my $bottom = $y2;

  my $left         = $xmid - $width/2;
  my $right        = $xmid + $width/2;

   my $fgcolor      = $self->fgcolor;

  my $span_height = $self->span_height;
  $self->filled_box($gd,$x1,$top,$x2,$top+$span_height);

  if ($self->top_spacing >= 8) {
    $top += $span_height;
    $gd->line($x1,$top+2,$x1,$top+4,$fgcolor);
    $gd->line($x2,$top+2,$x2,$top+4,$fgcolor);
    $gd->line($x1,$top+4,$left,$y1-4,$fgcolor);
    $gd->line($x2,$top+4,$right,$y1-4,$fgcolor);
    $gd->line($left,$y1-4,$left,$y1-2,$fgcolor);
    $gd->line($right,$y1-4,$right,$y1-2,$fgcolor);
  }

  $self->draw_contents($gd,$left,$y1,$right,$y2);

  my $pl = $self->pad_left;
  $self->draw_label($gd,$dx-$pl,$dy)       if $self->option('label');
  $self->draw_description($gd,$dx-$pl,$dy) if $self->option('description');
}

sub draw_contents {
  my $self = shift;
  my ($gd,$left,$top,$right,$bottom) = @_;
  $self->filled_box($gd,$left,$top,$right,$bottom);
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::fixedwidth - A base class fixed width glyphs

=head1 SYNOPSIS

 use Bio::Graphics;
 use Bio::Seq;
 use Bio::SeqFeature::Generic;

 my $bsg = 'Bio::SeqFeature::Generic';

 my $seq    = Bio::Seq->new(-length=>1000);

 my $whole  = $bsg->new(-display_name => 'Clone82',
 		        -start        => 1,
		        -end          => $seq->length);

 my $f1     = $bsg->new(-start        => 100,
		        -end          => 300,
		        -display_name => 'feature 1',
		       );

 my $f2      = $bsg->new(-start        => 500,
		         -end          => 800,
		         -display_name => 'feature 2',
		       );

 my $panel = Bio::Graphics::Panel->new(-length    => $seq->length,
				       -width     => 800,
				       -key_style => 'between',
				       -pad_left  => 10,
				       -pad_right => 10,
				      );

 $panel->add_track($whole,
		   -glyph    => 'arrow',
		   -double   => 1,
		   -tick     => 2,
		   -label    => 1,
		   );

 $panel->add_track([$f1,$f2],
		   -glyph    => 'fixedwidth',
		   -label    => 1,
                   -fixed_height => 20,
                   -fixed_width  => 20,
		   -key       => 'fixed width');

 binmode STDOUT;
 print $panel->png;

=head1 DESCRIPTION

This glyph is a base class for glyphs that wish to draw a fixed width
content, such as an icon, image, scatterplot, and it would be
inappropriate for the content to be stretched to match the start and
end point of the associated feature. Instead the glyph draws a simple
box spanning the feature's start:end region, two diagonal connecting
lines, and then a fixed width rectangle beneath the box.

This glyph does nothing very interesting itself. It is intended that
subclasses should override the draw_contents() method to draw
something interesting. See "Subclassing" for a simple example.

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

The following additional options are available to the "fixedwidth" glyph:

  Option            Description                       Default
  ------            -----------                       -------

  -fixed_width      Width of the content                 0


  -fixed_height     Height of the content                Same as -height

  -fixed_gap        Vertical gap between the box         8
                    that shows the extent of the
                    feature and the fixed-width
                    content.

                    If -fixed_gap is less than 8
                    then the diagonal connecting
                    lines are not drawn.

=head2 EXAMPLE SUBCLASS

To draw something interesting in the fixed rectangle, override the
draw_contents method. It takes four arguments consisting of the GD
object, and the left, top, right and bottom coordinates of the fixed
rectangle. Example:

 package Bio::Graphics::Glyph::fixedexample;
 use strict;
 use base 'Bio::Graphics::Glyph::fixedwidth';

 sub draw_contents {
   my $self = shift;
   my ($gd,$left,$top,$right,$bottom) = @_;
   $self->unfilled_box($gd,$left,$top,$right,$bottom);
   $gd->line($left,$top,$right,$bottom,$self->fgcolor);
   $gd->line($left,$bottom,$right,$top,$self->fgcolor);
 }

 1;

This will draw the outline of the fixed rectangle. The rectangle will
contain two diagonal lines. Not very interesting, but an example,
nonetheless.

See the stackedplot glyph for a more interesting subclass.

=head1 BUGS AND LIMITATIONS

This glyph should used as the base for the image glyph, but
isn't. This will be fixed.

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

Copyright (c) 2007 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
