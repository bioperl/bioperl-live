package Bio::Graphics::Glyph::spectrogram;
# $Id: spectrogram.pm,v 1.7 2008/09/20 15:46:39 lstein Exp $

use strict;
use Bio::Graphics::Glyph::generic;
use GD::Simple;
use GD;
use List::Util qw/sum max/;

use Data::Dumper;

use vars '@ISA';
@ISA = 'Bio::Graphics::Glyph::generic';

# saturation value is fixed at the maximum
use constant SAT  => 255;
use constant GDS  => GD::Simple->new(1,1);

# each spectrogram feature will be a standalone
# (unaggregated) feature.
sub draw {
  my $self = shift;
  my $gd = shift;
  my ( $x1, $y1, $x2, $y2 ) = $self->bounds(@_);
  my $v_offset = $y1;
  my $last_y   = 0;
  my $last_label;
  my $height   = $self->option('height');
  my $win      = $self->option('win');
  my $feat     = $self->feature;

  # API change?
  my $att;
  if (ref $feat->attributes) {
    $att = $feat->attributes;
  } 
  else {
    $att = {$feat->attributes};
  }

  my $max      = $att->{max};
  my $black    = $gd->colorResolve(0,0,0);

  my %seen;

  my $rows = @{$att->{g}};
  my $step = int $height/$rows;
  my $lbl_type = shift @{$att->{labels}} if defined $att->{labels};

  for my $g (@{$att->{g}}) {
    my $a   = shift @{$att->{a}};
    my $t   = shift @{$att->{t}};
    my $c   = shift @{$att->{c}};
    my $lbl = shift @{$att->{labels}} if defined $att->{labels}; 

    my ($hue,$bri) = get_bg_color($g,$a,$t,$c);
    $hue = int(($hue/360)  * 255);
    $hue += 255 if $hue < 0;
    $hue -= 255 if $hue > 255;
    $bri = int(($bri/$max) * 255);
    my @rgb = GDS->HSVtoRGB($hue,SAT,$bri);
    my $bgcolor = $gd->colorResolve(@rgb);
    $self->filled_box($gd, $x1, $y1, $x2, $y1+$step, $bgcolor, $bgcolor);

    if ($lbl && $y1 > $last_y+10) {

      my $label = sprintf '%4s', int $lbl;

      # print labels for the number closest to an integer
      # we use the previous value to catch the transition
      if ($last_label ne $label) {
	  $gd->string(gdSmallFont, $x1-25, $last_y-5, $last_label, $black);
      }
      # this will hopefully catch the label at the bottom of the stack
      elsif (!$att->{labels}->[0]) {
	$gd->string(gdSmallFont, $x1-25, $y1-5, $label, $black);
      }
      
      $last_y = $y1 + 15;
      $last_label = $label;
    }
  
    if ($lbl) {
      $gd->stringUp(gdSmallFont, $x1-27, $y2-5, $lbl_type, $black);
    }

    $y1 += $step;
  }
}


# HSV color space:
# Hue        (0-360 degrees)
# Saturation (0-100)
# Brightness (0-100)

sub get_bg_color {
  my ($g,$a,$t,$c) = @_;
  my $total = sum(@_) || return (0,0);
  my $max = max (@_);
  my $angle;

  # Assign the angular coordinate for the 
  # dominant base (>50% of signal) 
  if ($max == $g && $max >= $total/2) {
    $angle = 60;  # yellow
  }
  elsif ($max == $a && $max >= $total/2) {
    $angle = 240; # blue
  }
  elsif ($max == $t &&  $max >= $total/2) {
    $angle = 0;   # red
  }
  elsif ($max == $c &&  $max >= $total/2) {
    $angle = 120; # green
  }

  # or else take the weighted average coordinate
  # This is not perfect, as the coordinates are not
  # equidistant, but most spots will fall into
  # the above category
  else {
    my $acg     = 60*$g + 120*$c + 240*$a;
    my $t_angle = $acg/($g + $c + $a) > 180 ? 360 : 1;
    $angle  = ($t_angle*$t + $acg)/$total;
  }

  $angle += 0.5;

  return (int($angle), $total);
}

# make sure bumping is off to get an aligned spectrogram
sub bump { 0 }

1;

=head1 NAME

Bio::Graphics::Glyph::spectrogram - The "spectrogram" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel>, L<Bio::Graphics::Glyph>
      and L<Bio::Graphics::Browser::Plugin::Spectrogram> 

=head1 DESCRIPTION

This glyph is designed to draw DNA spectrograms for the
Spectrogram plugin.  It is not meant to be used as a
standalone glyph and has few public options.  Most of the
glyph's behavior is controlled via the spectrogram plugin.  

The glyph expects unaggregated 1D spectrogram features, each of which
is a vertical column, with one row for each integer frequency.
The number of frequencies is controlled by the window size and/or
the Spectrogram plugin.  The width of the feature corresponds to
the size of the overlap between adjacent windows. 
The values for each frequency are in four channels,
one for each base.  The color of each row ("spot")
represents the dominant base(s) and the intensity represents 
the magnitude of the signal at that frequency

The entire 2D spectrogram is a series of 
unaggregated, unbumped 1D spectrogram features.
 
The spectrogram glyph assigns colors using the HSV color space,
where an angular coordinate for hue is assigned to each base 
(G yellow [60]; A blue [240]; C green [120]; T red[0/360]). 
The saturation value is fixed at the maximum of 100 and the brightness 
value is scaled according to the magnitude for each frequency,
ranging from black to the pure hue.

The hue is determined in one of two ways:  

If the signal for one base is dominant (> 50% of total for the four
channels) the angular coordinate for that base is used.
The brightness is calculated using the total signal from
all four channels.

If no base has a dominant signal, the weighted average angular
coordinate is calculated using the relative contribution
from each channel.  The brightness is calculated from
the total signal from all four channels.

The y-axis labels require at least 40 pixels of left-padding.
They will be truncated if less than 40 of padding is specified
in the configuration file.


=head2 OPTIONS

The following standard options are accepted:

  Option      Description                      Default
  ------      -----------                      -------

  -height     Height of glyph		       calculated

  -bump       Whether to bump features         off
 
The following glyph-specific options are also used:
  
  -win        window size used to calculate    calculated
              the spectrogram values


=head1 BUGS

Please report them.

=head1 AUTHOR

Sheldon McKay E<lt>mckays@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut
