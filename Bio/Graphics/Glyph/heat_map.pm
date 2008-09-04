package Bio::Graphics::Glyph::heat_map;
#$Id: heat_map.pm,v 1.4.2.3 2007/10/17 01:48:22 lstein Exp $

use strict;
use Bio::Graphics::Glyph::minmax;

# A glyph to draw a heat map for scored features along a continuous color
# gradient calculated in HSV color space

use vars '@ISA';
@ISA = qw/Bio::Graphics::Glyph::minmax/;

# set up getter/setter methods
BEGIN {
  no strict 'refs';

  my @subs = qw/ h_start   s_start   v_start h_range s_range  v_range
                 min_score max_score low_rgb low_hsv high_rgb score_range/;

  for my $sub ( @subs ) {
    *{$sub} = sub {
      my ($self, $v) = @_;
      my $k = "_$sub";

      if (defined $v) {
	$self->{$k} = $v;
      }

      return $self->{$k};
    }
  }
}

sub draw {
  my $self = shift;

  my @parts = $self->parts;
  @parts    = $self if !@parts && $self->level == 0;
  return $self->SUPER::draw(@_) unless @parts;

  $self->calculate_gradient(\@parts);
  my $low_rgb = $self->low_rgb;

  for my $part (@parts) {
    my $score = $part->feature->score;

    # use start color if no score or no score gradient
    unless (defined $score && $self->score_range ) {
      $part->{partcolor} = $self->color_index(@$low_rgb);
    }
    else {
      my @rgb = $self->calculate_color($score);
      $part->{partcolor} = $self->color_index(@rgb);
    }
    
  }

  return $self->SUPER::draw(@_);
}

# We want an exact match, so allocate the color
# if required
sub color_index {
  my ($self, @rgb) = @_;
  my $gd = $self->panel->gd;
  return $gd->colorResolve(@rgb);
}

# Override minmax method to get user supplied
# values.  This will be helpful for single or
# unaggregated features.
sub minmax {
  my ($self, $parts) = @_;
  my $min  = $self->option('min_score');
  my $max  = $self->option('max_score');
  return ($min,$max) if $min && $max && $min < $max;
  return (0,$max)    if $max && !$min; # minscore may be zero
  return (0,100)     unless $parts;
  return $self->SUPER::minmax($parts);
}

# convert named color or hex string to RGB value, then HSV 
sub color2hsv {
  my ($self,$color) = @_;
  my $color_idx = $self->panel->translate_color($color);
  my @rgb = $self->panel->rgb($color_idx);
  return [$self->RGBtoHSV(@rgb)];
}

sub calculate_gradient {
  my ($self, $parts) = @_;
  my $start_color = lc $self->option('start_color') || 'white';
  my $stop_color  = lc $self->option('end_color')   || 'red';
  my $hsv_start   = $self->color2hsv($start_color);
  my $hsv_stop    = $self->color2hsv($stop_color);

  my ($h_start,$s_start,$v_start) = @$hsv_start;
  my ($h_stop,$s_stop,$v_stop )   = @$hsv_stop;
  my $h_range = $h_stop - $h_start;
  my $s_range = $s_stop - $s_start;
  my $v_range = $v_stop - $v_start;

  # override brightness and saturation if required
  if (my $bri = $self->option('brightness')) {
    $bri = int($bri*255/100 + 0.5);
    $v_start = $v_stop = $bri;
    $v_range = 0;
  }
  if (my $sat = $self->option('saturation')) {
    $sat = int($sat*255/100 + 0.5);
    $s_start = $s_stop = $sat;
    $s_range = 0;
  }
  if ($self->option('pure_hue')) {
    $hsv_start = [$h_start,255,255];
    $hsv_stop  = [$h_stop,255,255];
    $v_start   = $v_stop  = 255;
    $s_start   = $s_stop  = 255;
    $v_range   = $s_range = 0;
  }

  # darkness or monochrome gradient?
  if ( !_isa_color($start_color) || !_isa_color($stop_color) ) {
    # hue (H) is fixed
    $h_range = 0;

    #    gradient         S       V    
    # white -> color    0->255   255
    # color -> white    255->0   255
    # white -> black    0        255->0
    # black -> white    0        0->255
    # black -> color    0->255   0->255
    # color -> black    255->0   255->0
    if ( $start_color eq 'white' && _isa_color($stop_color) ) {
      $s_range = 255;
      $s_start = 0;
      $v_range = 0;
      $v_start = 255;
      $h_start = $h_stop;
    }
    elsif ( _isa_color($start_color) && $stop_color eq 'white' ) {
      $s_range = -255;
      $s_start = 255;
      $v_range = 0;
      $v_start = 255;
    }
    elsif ( $start_color eq 'white' ) { # end black
      $s_range = 0;
      $s_start = 0;
      $v_range = -255;
      $v_start = 255;
    }
    elsif ( $stop_color eq 'white' ) { # start black
      $s_range = 0;
      $s_start = 0;
      $v_range = 255;
      $v_start = 0;
    }
    elsif ( _isa_color($start_color) ) { # end black
      $s_range = 255;
      $s_start = 0;
      $v_range = 255;
      $v_start = 0;
    }
    elsif ( _isa_color($stop_color) ) { # start black
      $s_range = -255;
      $s_start = 255;
      $v_range = -255;
      $v_start = 255;
    }
	
  }

  # store gradient info
  $self->h_range($h_range);
  $self->h_start($h_start);
  $self->s_start($s_start);
  $self->v_start($v_start);
  $self->s_range($s_range);
  $self->v_range($v_range);

  # store score info
  my ($min,$max) = $self->minmax($parts);
  $self->score_range($max - $min);
  $self->min_score($min);
  $self->max_score($max);
  
  # store color extremes
  my @low_rgb  = $self->HSVtoRGB(@$hsv_start);
  my @high_rgb = $self->HSVtoRGB(@$hsv_stop);
  $self->low_hsv($hsv_start);
  $self->high_rgb(\@high_rgb);
  $self->low_rgb(\@low_rgb);
  return 1;
}

sub _isa_color {
  my $color = shift;
  return $color =~ /white|black|FFFFFF|000000/i ? 0 : 1;
}

sub calculate_color {
  my ($self,$score) = @_;
  $score ||= 0;

  # relative score
  my $min   = $self->min_score;
  my $max   = $self->max_score;
  my $range = $self->score_range;

  # reset off-scale scores
  $score = $min if $score < $min;
  $score = $max if $score > $max;
  my $score_diff = ($score - $min)/$range;

  # Hue 
  my $hue    = $self->h_start;
  my $h_diff = $score_diff * $self->h_range;
  $hue += $h_diff;
  $hue = int($hue+0.5);

  # Saturation
  my $sat = $self->s_start;
  $sat += $score_diff * $self->s_range; 
  $sat = int($sat+0.5);

  # Brightness
  my $bri = $self->v_start;
  $bri += $score_diff * $self->v_range;
  $bri = int($bri + 0.5);

  return $self->HSVtoRGB($hue,$sat,$bri);
}

# synthesize a key glyph
sub keyglyph {
  my $self = shift;
  my $scale = 1/$self->scale;  # base pairs/pixel
  my $offset = $self->panel->offset;
  my ($min,$max) = $self->minmax;
  my $range = $max - $min;
  my ($segments, $low);

  for my $start (0..9) {
    $start *= 10;
    push @$segments, [ $start*$scale + $offset, ($start + 10)*$scale + $offset ];
  }

  my $feature = Bio::Graphics::Feature->new( -segments => $segments,
					     -name     => $self->option('key'),
					     -strand   => '+1' );

  for (0..9) {
    my $score += ($range/10) * $_; 
    ($feature->segments)[$_]->score($score);
  }

  my $factory = $self->factory->clone;
  $factory->set_option(label => 1);
  $factory->set_option(bump  => 0);
  $factory->set_option(min_score  => 0);
  $factory->set_option(max_score  => 100);
  return $factory->make_glyph(0,$feature);
}

sub bgcolor { 
  my $self = shift;
  return defined $self->{partcolor} ? $self->{partcolor} : $self->SUPER::bgcolor;
}
sub fgcolor {
  my $self = shift;
  return $self->bgcolor;
}

sub RGBtoHSV {
  my ($self, $r, $g ,$bl) = @_;
  my ($min,undef,$max) = sort {$a<=>$b} ($r,$g,$bl);

  my $range = $max - $min or return (0,0,$r);
  my $v = $max;
  my $s = 255 * ($max - $min)/$max;
  my $h;
  
  if ($max == $r) {
    $h = 60 * ($g-$bl)/$range;
  }
  elsif ($max == $g) {
    $h = 60 * ($bl-$r)/$range + 120;
  }
  else {
    $h = 60 * ($r-$g)/$range + 240;
  }

  $h  = int($h*255/360 + 0.5);
  $h += 255 if $h < 0;
  $h -= 255 if $h > 255;

  return ($h, $s, $v);
}

# method courtesy of Lincoln Stein
sub HSVtoRGB {
  my $self = shift;
  @_ == 3 or die "Usage: GD::Simple->HSVtoRGB(\$hue,\$saturation,\$value)";

  my ($h,$s,$v)=@_;
  my ($r,$g,$b,$i,$f,$p,$q,$t);

  if( $s == 0 ) {
    ## achromatic (grey)
    return ($v,$v,$v);
  }
  $h %= 255;
  $s /= 255;                      ## scale saturation from 0.0-1.0
  $h /= 255;                      ## scale hue from 0 to 1.0
  $h *= 360;                      ## and now scale it to 0 to 360

  $h /= 60;                       ## sector 0 to 5
  $i = $h % 6;
  $f = $h - $i;                   ## factorial part of h
  $p = $v * ( 1 - $s );
  $q = $v * ( 1 - $s * $f );
  $t = $v * ( 1 - $s * ( 1 - $f ) );

  if($i<1) {
    $r = $v;
    $g = $t;
    $b = $p;
  } elsif($i<2){
    $r = $q;
    $g = $v;
    $b = $p;
  } elsif($i<3){
    $r = $p;
    $g = $v;
    $b = $t;
  } elsif($i<4){
    $r = $p;
    $g = $q;
    $b = $v;
  } elsif($i<5){
    $r = $t;
    $g = $p;
    $b = $v;
  } else {
    $r = $v;
    $g = $p;
    $b = $q;
  }
  return (int($r+0.5),int($g+0.5),int($b+0.5));
}


1;

=head1 NAME

Bio::Graphics::Glyph::heat_map - The "heat_map" glyph

=head1 SYNOPSIS

See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws "scored" features using a continuous
color gradient is the HSV color space. The color of 
each segment is proportionate to the score.

=head1 OPTIONS

=head2 Global glyph options:

See L<Bio::Graphics::Glyph::generic>

=head2 Glyph-specific options:

The color_map glyph recognizes the following
glyph-specific options:

  Option      Description                   Default
  ------      -----------                   -------

  -start_color Beginning of the color       white
               gradient, expressed as a 
               named color or RGB hex 
               string
 
  -end_color   End of the color gradient    red

  -brightness  Color brilliance:  0-100     Calculated
               This will override the 
               value from the named
               color

  -saturation  Color saturation: 0-100      Calculated
               This will override the
               value from the named
               color

  -pure_hue    Use the pure hue (bright-    0 (false)
               ness and saturation both
               at 100) for the named coloe
               
  -max_score   Maximum value of the	    Calculated
               feature's "score" attribute

  -min_score   Minimum value of the         Calculated
               feature's "score" attribute

  -vary_fg     Vary the foreground color    1 (true)
               with the background color



If max_score and min_score are not specified, the glyph will
calculate the local maximum and minimum scores at run time.  If single
features, unaggregated features, or multiple aggregates are being drawn, 
this will result in an inconsistent color scale.  It is recommended
that global maximum and minimum scores be specified in the track 
configuration. Since many scoring functions are exponential,
you may wish to take the log of your scores before passing
them to this glyph.

=head2 Color Gradients

The color gradient is calculating by
progressing along the rainbow spectrum from red to violet,
also incrementing brightness and saturation, all proportate 
to the score value.  To vary the hue only, "pure" hues can
be used.  Pure hues have brightness and saturation values of
100. Some examples, in order, are red, yellow, lime, 
aqua/cyan, blue and magenta.  The gradient can progress in
reverse orientation with the respect to the visible light
spectrum if a lower-order color is used as the start and a higher
order color used as the end (for example lime->red).

Using the "pure_hue" option results in a brighter, more vibrant
color spectrum, Choosing darker start and end colors, such as
green or maroon, will result in a darker spectrum.  
A single color spectrum can be created by using black
or white as the start or end color.  

A grayscale spectrum will result if black and white 
are used as start and end colors.  One example of an
effective visual heat map is to progress from 
white->red.

For the start_color and end_color options, 140 named webcolors
and their corresponsing RGB hex codes (listed below) are supported.

 steelblue           	#4682B4
 royalblue           	#041690
 cornflowerblue      	#6495ED
 lightsteelblue      	#B0C4DE
 mediumslateblue     	#7B68EE
 slateblue           	#6A5ACD
 darkslateblue       	#483D8B
 midnightblue        	#191970
 navy                	#000080
 darkblue            	#00008B
 mediumblue          	#0000CD
 blue                	#0000FF
 dodgerblue          	#1E90FF
 deepskyblue         	#00BFFF
 lightskyblue        	#87CEFA
 skyblue             	#87CEEB
 lightblue           	#ADD8E6
 powderblue          	#B0E0E6
 azure               	#F0FFFF
 lightcyan           	#E0FFFF
 paleturquoise       	#AFEEEE
 mediumturquoise     	#48D1CC
 lightseagreen       	#20B2AA
 darkcyan            	#008B8B
 teal                	#008080
 cadetblue           	#5F9EA0
 darkturquoise       	#00CED1
 aqua                	#00FFFF
 cyan                	#00FFFF
 turquoise           	#40E0D0
 aquamarine          	#7FFFD4
 mediumaquamarine    	#66CDAA
 darkseagreen        	#8FBC8F
 mediumseagreen      	#3CB371
 seagreen            	#2E8B57
 darkgreen           	#006400
 green               	#008000
 forestgreen         	#228B22
 limegreen           	#32CD32
 lime                	#00FF00
 chartreuse          	#7FFF00
 lawngreen           	#7CFC00
 greenyellow         	#ADFF2F
 yellowgreen         	#9ACD32
 palegreen           	#98FB98
 lightgreen          	#90EE90
 springgreen         	#00FF7F
 mediumspringgreen   	#00FA9A
 darkolivegreen      	#556B2F
 olivedrab           	#6B8E23
 olive               	#808000
 darkkhaki           	#BDB76B
 darkgoldenrod       	#B8860B
 goldenrod           	#DAA520
 gold                	#FFD700
 yellow              	#FFFF00
 khaki               	#F0E68C
 palegoldenrod       	#EEE8AA
 blanchedalmond      	#FFEBCD
 moccasin            	#FFE4B5
 wheat               	#F5DEB3
 navajowhite         	#FFDEAD
 burlywood           	#DEB887
 tan                 	#D2B48C
 rosybrown           	#BC8F8F
 sienna              	#A0522D
 saddlebrown         	#8B4513
 chocolate           	#D2691E
 peru                	#CD853F
 sandybrown          	#F4A460
 darkred             	#8B0000
 maroon              	#800000
 brown               	#A52A2A
 firebrick           	#B22222
 indianred           	#CD5C5C
 lightcoral          	#F08080
 salmon              	#FA8072
 darksalmon          	#E9967A
 lightsalmon         	#FFA07A
 coral               	#FF7F50
 tomato              	#FF6347
 darkorange          	#FF8C00 
 orange              	#FFA500
 orangered           	#FF4500
 crimson             	#DC143C
 red                 	#FF0000
 deeppink            	#FF1493
 fuchsia             	#FF00FF
 magenta             	#FF00FF
 hotpink             	#FF69B4
 lightpink           	#FFB6C1
 pink                	#FFC0CB
 palevioletred       	#DB7093
 mediumvioletred     	#C71585
 purple              	#800080
 darkmagenta         	#8B008B
 mediumpurple        	#9370DB
 blueviolet          	#8A2BE2
 indigo              	#4B0082
 darkviolet          	#9400D3
 darkorchid          	#9932CC
 mediumorchid        	#BA55D3 
 orchid              	#DA70D6 
 violet              	#EE82EE
 plum                	#DDA0DD
 thistle             	#D8BFD8
 lavender            	#E6E6FA
 ghostwhite          	#F8F8FF
 aliceblue           	#F0F8FF
 mintcream           	#F5FFFA
 honeydew            	#F0FFF0
 lightgoldenrodyellow	#FAFAD2
 lemonchiffon        	#FFFACD
 cornsilk            	#FFF8DC
 lightyellow         	#FFFFE0
 ivory               	#FFFFF0
 floralwhite         	#FFFAF0
 linen               	#FAF0E6
 oldlace             	#FDF5E6
 antiquewhite        	#FAEBD7
 bisque              	#FFE4C4
 peachpuff           	#FFDAB9
 papayawhip          	#FFEFD5
 beige               	#F5F5DC
 seashell            	#FFF5EE
 lavenderblush       	#FFF0F5
 mistyrose           	#FFE4E1
 snow                	#FFFAFA
 white               	#FFFFFF
 whitesmoke          	#F5F5F5
 gainsboro           	#DCDCDC
 lightgrey           	#D3D3D3
 silver              	#C0C0C0
 darkgray            	#A9A9A9
 gray                	#808080
 lightslategray      	#778899
 slategray           	#708090
 dimgray             	#696969
 darkslategray       	#2F4F4F
 black               	#000000


=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Glyph>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::heterogeneous_segments>,
L<Bio::Graphics::Feature>,
L<Bio::DB::GFF>
L<GD>        

=head1 AUTHOR

Sheldon McKay E<lt>mckays@cshl.eduE<gt>

Copyright (c) 2006 Cold Spring Harbor Laboratory

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut
