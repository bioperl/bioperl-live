package Bio::Graphics::Panel;

use strict;
use Bio::Graphics::Glyph::Factory;
use Bio::Graphics::Feature;
use Carp 'cluck';
use GD;
use vars '$VERSION';

$VERSION = '1.01';

use constant KEYLABELFONT => gdMediumBoldFont;
use constant KEYSPACING   => 5; # extra space between key columns
use constant KEYPADTOP    => 5;  # extra padding before the key starts
use constant KEYCOLOR     => 'wheat';
use constant KEYSTYLE     => 'bottom';
use constant KEYALIGN     => 'left';
use constant GRIDCOLOR    => 'lightcyan';

my %COLORS;  # translation table for symbolic color names to RGB triple

# Create a new panel of a given width and height, and add lists of features
# one by one
sub new {
  my $class = shift;
  my %options = @_;

  $class->read_colors() unless %COLORS;

  my $length = $options{-length} || 0;
  my $offset = $options{-offset}  || 0;
  my $spacing = $options{-spacing} || 5;
  my $bgcolor = $options{-bgcolor} || 0;
  my $keyfont = $options{-key_font} || KEYLABELFONT;
  my $keycolor = $options{-key_color} || KEYCOLOR;
  my $keyspacing = $options{-key_spacing} || KEYSPACING;
  my $keystyle = $options{-key_style} || KEYSTYLE;
  my $keyalign = $options{-key_align} || KEYALIGN;
  my $allcallbacks = $options{-all_callbacks} || 0;
  my $gridcolor    = $options{-gridcolor} || GRIDCOLOR;
  my $grid         = $options{-grid} || 0;

  $offset   ||= $options{-segment}->start-1 if $options{-segment};
  $length   ||= $options{-segment}->length  if $options{-segment};

  $offset   ||= $options{-start}-1 if defined $options{-start};
  $length   ||= $options{-stop}-$options{-start}+1 
     if defined $options{-start} && defined $options{-stop};

  return bless {
		tracks => [],
		width      => $options{-width} || 600,
		pad_top    => $options{-pad_top}||0,
		pad_bottom => $options{-pad_bottom}||0,
		pad_left   => $options{-pad_left}||0,
		pad_right  => $options{-pad_right}||0,
		length => $length,
		offset => $offset,
		gridcolor => $gridcolor,
		grid    => $grid,
		bgcolor => $bgcolor,
		height => 0, # AUTO
		spacing => $spacing,
		key_font => $keyfont,
		key_color => $keycolor,
		key_spacing => $keyspacing,
		key_style => $keystyle,
		key_align => $keyalign,
		all_callbacks => $allcallbacks,
	       },$class;
}

sub pad_left {
  my $self = shift;
  my $g = $self->{pad_left};
  $self->{pad_left} = shift if @_;
  $g;
}
sub pad_right {
  my $self = shift;
  my $g = $self->{pad_right};
  $self->{pad_right} = shift if @_;
  $g;
}
sub pad_top {
  my $self = shift;
  my $g = $self->{pad_top};
  $self->{pad_top} = shift if @_;
  $g;
}
sub pad_bottom {
  my $self = shift;
  my $g = $self->{pad_bottom};
  $self->{pad_bottom} = shift if @_;
  $g;
}

# numerous direct calls into array used here for performance considerations
sub map_pt {
  my $self   = shift;
  my $offset = $self->{offset};
  my $scale  = $self->{scale} || $self->scale;
  my $pl = $self->{pad_left};
  my $pr = $self->{width} - $self->{pad_right};
  my @result;
  foreach (@_) {
    my $val = int (0.5 + $pl + ($_-$offset-1) * $scale);
    $val = $pl-1 if $val < $pl;
    $val = $pr+1 if $val > $pr;
    push @result,$val;
  }
  @result;
}

sub map_no_trunc {
  my $self   = shift;
  my $offset = $self->{offset};
  my $scale  = $self->{scale} || $self->scale;
  my $pl = $self->{pad_left};
  my $pr = $self->{width} - $self->{pad_right};
  my @result;
  foreach (@_) {
    my $val = int (0.5 + $pl + ($_-$offset-1) * $scale);
    push @result,$val;
  }
  @result;
}

sub scale {
  my $self = shift;
  $self->{scale} ||= ($self->{width}-$self->pad_left-$self->pad_right-1)/($self->length-1);
}

sub start { shift->{offset}+1}
sub end   { $_[0]->start + $_[0]->{length}-1}

sub offset { shift->{offset} }
sub width {
  my $self = shift;
  my $d = $self->{width};
  $self->{width} = shift if @_;
  $d;
#  $d + $self->pad_left + $self->pad_right;
}

sub left {
  my $self = shift;
  $self->pad_left;
}
sub right {
  my $self = shift;
  $self->width - $self->pad_left;
}

sub spacing {
  my $self = shift;
  my $d = $self->{spacing};
  $self->{spacing} = shift if @_;
  $d;
}

sub key_spacing {
  my $self = shift;
  my $d = $self->{key_spacing};
  $self->{key_spacing} = shift if @_;
  $d;
}

sub length {
  my $self = shift;
  my $d = $self->{length};
  if (@_) {
    my $l = shift;
    $l = $l->length if ref($l) && $l->can('length');
    $self->{length} = $l;
  }
  $d;
}

sub gridcolor {shift->{gridcolor}}

sub all_callbacks { shift->{all_callbacks} }

sub add_track {
  my $self = shift;
  $self->_do_add_track(scalar(@{$self->{tracks}}),@_);
}

sub unshift_track {
  my $self = shift;
  $self->_do_add_track(0,@_);
}

sub insert_track {
  my $self = shift;
  my $position = shift;
  $self->_do_add_track($position,@_);
}


# create a feature and factory pair
# see Factory.pm for the format of the options
# The thing returned is actually a generic Glyph
sub _do_add_track {
  my $self     = shift;
  my $position = shift;

  # due to indecision, we accept features
  # and/or glyph types in the first two arguments
  my ($features,$glyph_name) = ([],undef);
  while ( @_ && $_[0] !~ /^-/) {
    my $arg = shift;
    $features   = $arg and next if ref($arg);
    $glyph_name = $arg and next unless ref($arg);
  }

  my %args = @_;
  my ($map,$ss,%options);

  foreach (keys %args) {
    (my $canonical = lc $_) =~ s/^-//;
    if ($canonical eq 'glyph') {
      $map = $args{$_};
      delete $args{$_};
    } elsif ($canonical eq 'stylesheet') {
      $ss  = $args{$_};
      delete $args{$_};
    } else {
      $options{$canonical} = $args{$_};
    }
  }

  $glyph_name = $map if defined $map;
  $glyph_name ||= 'generic';

  my $panel_map = ref($map) eq 'CODE'
    ?  sub {
          my $feature = shift;
	  return 'track' if eval { $feature->primary_tag  eq 'track' };
	  return 'group' if eval { $feature->primary_tag  eq 'group' };
	  return $map->($feature);
	}
      :
	sub {
	  my $feature = shift;
	  return 'track' if eval { $feature->primary_tag  eq 'track' };
	  return 'group' if eval { $feature->primary_tag  eq 'group' };
	  return $glyph_name;
	};

  $self->_add_track($position,$features,-map=>$panel_map,-stylesheet=>$ss,-options=>\%options);
}

sub _add_track {
  my $self = shift;
  my ($position,$features,@options) = @_;

  # build the list of features into a Bio::Graphics::Feature object
  $features = [$features] unless ref $features eq 'ARRAY';

  # optional middle-level glyph is the group
  foreach my $f (grep {ref $_ eq 'ARRAY'} @$features) {
    next unless ref $f eq 'ARRAY';
    $f = Bio::Graphics::Feature->new(
				     -segments=>$f,
				     -type => 'group'
				    );
  }

  # top-level glyph is the track
  my $feature = Bio::Graphics::Feature->new(
					    -segments=>$features,
					    -start => $self->offset+1,
					    -stop  => $self->offset+$self->length,
					    -type => 'track'
					   );

  my $factory = Bio::Graphics::Glyph::Factory->new($self,@options);
  my $track   = $factory->make_glyph(-1,$feature);

  splice(@{$self->{tracks}},$position,0,$track);
  return $track;
}

sub height {
  my $self = shift;
  my $spacing    = $self->spacing;
  my $key_height = $self->format_key;
  my $height = 0;
  foreach (@{$self->{tracks}}) {
    next unless $_->parts;
    $height += $_->layout_height + $spacing;
  }
  return $height + $key_height + $self->pad_top + $self->pad_bottom;
}

sub gd {
  my $self = shift;

  return $self->{gd} if $self->{gd};

  my $width  = $self->width;
  my $height = $self->height;

  my $gd = GD::Image->new($width,$height);
  my %translation_table;
  for my $name ('white','black',keys %COLORS) {
    my $idx = $gd->colorAllocate(@{$COLORS{$name}});
    $translation_table{$name} = $idx;
  }

  $self->{translations} = \%translation_table;
  $self->{gd}                = $gd;
  $gd->fill(0,0,$self->bgcolor) if $self->bgcolor;

  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $offset = $pt;
  my $keyheight   = $self->{key_font}->height;
  my $between_key = $self->{key_style} eq 'between';
  my $spacing = $self->spacing;

  # we draw in two steps, once for background of tracks, and once for
  # the contents.  This allows the grid to sit on top of the track background.
  for my $track (@{$self->{tracks}}) {
    next unless $track->parts;
    $gd->filledRectangle($pl,
			 $offset,
			 $width-$self->pad_right,
			 $offset+$track->layout_height
			 + ($self->{key_style} eq 'between' ? $self->{key_font}->height : 0),
			 $track->tkcolor)
      if defined $track->tkcolor;
    $offset += $keyheight if $between_key && $track->option('key');
    $offset += $track->layout_height + $spacing;
  }

  $self->draw_grid($gd)  if $self->{grid};

  $offset = $pt;
  for my $track (@{$self->{tracks}}) {
    next unless $track->parts;
    $offset += $self->draw_between_key($gd,$track,$offset) if $between_key && $track->option('key');
    $track->draw($gd,0,$offset,0,1);
    $offset += $track->layout_height + $spacing;
  }

  $self->draw_bottom_key($gd,$pl,$offset) if $self->{key_style} eq 'bottom';
  return $self->{gd} = $gd;
}

sub boxes {
  my $self = shift;
  my @boxes;
  my $offset = 0;
  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $between = $self->{key_style} eq 'between';
  for my $track (@{$self->{tracks}}) {
    next unless $track->parts;
    $offset += $self->{key_font}->height if $between && $track->option('key');
    my $boxes = $track->boxes(0,$offset+$pt);
    push @boxes,@$boxes;
    $offset += $track->layout_height + $self->spacing;
  }
  return wantarray ? @boxes : \@boxes;
}

# draw the keys -- between
sub draw_between_key {
  my $self   = shift;
  my ($gd,$track,$offset) = @_;
  my $key = $track->option('key') or return 0;
  my $x =   $self->{key_align} eq 'center' ? $self->width - (CORE::length($key) * $self->{key_font}->width)/2
          : $self->{key_align} eq 'right'  ? $self->width - CORE::length($key)
          : $self->pad_left;
  $gd->string($self->{key_font},$x,$offset,$key,1);
  return $self->{key_font}->height;
}

# draw the keys -- bottom
sub draw_bottom_key {
  my $self = shift;
  my ($gd,$left,$top) = @_;
  my $key_glyphs = $self->{key_glyphs} or return;

  my $color = $self->translate_color($self->{key_color});
  $gd->filledRectangle($left,$top,$self->width - $self->pad_right,$self->height-$self->pad_bottom,$color);
  $gd->string($self->{key_font},$left,KEYPADTOP+$top,"KEY:",1);
  $top += $self->{key_font}->height + KEYPADTOP;

  $_->draw($gd,$left,$top) foreach @$key_glyphs;
}

# Format the key section, and return its height
sub format_key {
  my $self = shift;

  return $self->{key_height} if defined $self->{key_height};

  if ($self->{key_style} eq 'between') {
    my @key_tracks = grep {$_->option('key')} @{$self->{tracks}};
    return $self->{key_height} = @key_tracks * $self->{key_font}->height;
  } elsif ($self->{key_style} eq 'bottom') {

    my ($height,$width) = (0,0);
    my %tracks;
    my @glyphs;

    # determine how many glyphs become part of the key
    # and their max size
    for my $track (@{$self->{tracks}}) {
      next unless $track->option('key');

      my $glyph;
      if (my @parts = $track->parts) {
	$glyph = $parts[0]->keyglyph;
      } else {
	my $t = Bio::Graphics::Feature->new(-segments=>
					    [Bio::Graphics::Feature->new(-start => $self->offset,
									 -stop  => $self->offset+$self->length)]);
	my $g = $track->factory->make_glyph(0,$t);
	$glyph = $g->keyglyph;
      }
      next unless $glyph;

      $tracks{$track} = $glyph;
      my ($h,$w) = ($glyph->layout_height,
		    $glyph->layout_width);
      $height = $h if $h > $height;
      $width  = $w if $w > $width;
      push @glyphs,$glyph;
    }

    $width += $self->key_spacing;

    # no key glyphs, no key
    return $self->{key_height} = 0 unless @glyphs;

    # now height and width hold the largest glyph, and $glyph_count
    # contains the number of glyphs.  We will format them into a
    # box that is roughly 3 height/4 width (golden mean)
    my $rows = 0;
    my $cols = 0;
    my $maxwidth = $self->width - $self->pad_left - $self->pad_right;
    while (++$rows) {
      $cols = @glyphs / $rows;
      $cols = int ($cols+1) if $cols =~ /\./;  # round upward for fractions
      my $total_width  = $cols * $width;
      my $total_height = $rows * $width;
      last if $total_width < $maxwidth;
    }

    # move glyphs into row-major format
    my $spacing = $self->key_spacing;
    my $i = 0;
    for (my $c = 0; $c < $cols; $c++) {
      for (my $r = 0; $r < $rows; $r++) {
	my $x = $c * ($width  + $spacing);
	my $y = $r * ($height + $spacing);
	next unless defined $glyphs[$i];
	$glyphs[$i]->move($x,$y);
	$i++;
      }
    }

    $self->{key_glyphs} = \@glyphs;     # remember our key glyphs
    # remember our key height
    return $self->{key_height} =
      ($height+$spacing) * $rows + $self->{key_font}->height +KEYPADTOP;
  }

  else {  # no known key style, neither "between" nor "bottom"
    return $self->{key_height} = 0;
  }
}

# draw a grid
sub draw_grid {
  my $self = shift;
  my $gd = shift;

  my $gridcolor = $self->translate_color($self->{gridcolor});
  my @positions;
  if (ref $self->{grid} eq 'ARRAY') {
    @positions = @{$self->{grid}};
  } else {
    my ($major,$minor) = $self->ticks;
    my $first_tick = $minor * int(0.5 + $self->start/$minor);
    for (my $i = $first_tick; $i < $self->end; $i += $minor) {
      push @positions,$i;
    }
  }
  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $pb = $self->height - $self->pad_bottom;
  for my $tick (@positions) {
    my ($pos) = $self->map_pt($tick);
    $gd->line($pos,$pt,$pos,$pb,$gridcolor);
  }
}

# calculate major and minor ticks, given a start position
sub ticks {
  my $self = shift;
  my ($length,$minwidth) = @_;

  $length   = $self->{length}       unless defined $length;
  $minwidth = gdSmallFont->width*7  unless defined $minwidth;

  my ($major,$minor);

  # figure out tick mark scale
  # we want no more than 1 major tick mark every 40 pixels
  # and enough room for the labels
  my $scale = $self->scale;

  my $interval = 1;

  while (1) {
    my $pixels = $interval * $scale;
    last if $pixels >= $minwidth;
    $interval *= 10;
  }
  return ($interval,$interval/10);
}

# reverse of translate(); given index, return rgb tripler
sub rgb {
  my $self = shift;
  my $idx  = shift;
  my $gd = $self->{gd} or return;
  return $gd->rgb($idx);
}

sub translate_color {
  my $self = shift;
  my @colors = @_;
  if (@colors == 3) {
    my $gd = $self->gd or return 1;
    return $gd->colorClosest(@colors);
  }
  elsif ($colors[0] =~ /^\#([0-9A-F]{2})([0-9A-F]{2})([0-9A-F]{2})$/i) {
    my $gd = $self->gd or return 1;
    my ($r,$g,$b) = (hex($1),hex($2),hex($3));
    return $gd->colorClosest($r,$g,$b);
  }
  else {
    my $color = $colors[0];
    my $table = $self->{translations} or return 1;
    return defined $table->{$color} ? $table->{$color} : 1;
  }
}

sub bgcolor {
   my $self = shift;
   return unless $self->{bgcolor};
   $self->translate_color($self->{bgcolor});
}

sub set_pen {
  my $self = shift;
  my ($linewidth,$color) = @_;
  return $self->{pens}{$linewidth,$color} if $self->{pens}{$linewidth,$color};

  my $pen = $self->{pens}{$linewidth} = GD::Image->new($linewidth,$linewidth);
  my @rgb = $self->rgb($color);
  my $bg = $pen->colorAllocate(255,255,255);
  my $fg = $pen->colorAllocate(@rgb);
  $pen->fill(0,0,$fg);
  $self->{gd}->setBrush($pen);
  return gdBrushed;
}

sub png {
  my $gd = shift->gd;
  $gd->png;
}

sub read_colors {
  my $class = shift;
  while (<DATA>) {
    chomp;
    last if /^__END__/;
    my ($name,$r,$g,$b) = split /\s+/;
    $COLORS{$name} = [hex $r,hex $g,hex $b];
  }
}

sub color_name_to_rgb {
  my $class = shift;
  my $color_name  = shift;
  $class->read_colors() unless %COLORS;
  return unless $COLORS{$color_name};
  return wantarray ? @{$COLORS{$color_name}}
                   : $COLORS{$color_name};
}

sub color_names {
    my $class = shift;
    $class->read_colors unless %COLORS;
    return wantarray ? keys %COLORS : [keys %COLORS];
}

1;

__DATA__
white                FF           FF            FF
black                00           00            00
aliceblue            F0           F8            FF
antiquewhite         FA           EB            D7
aqua                 00           FF            FF
aquamarine           7F           FF            D4
azure                F0           FF            FF
beige                F5           F5            DC
bisque               FF           E4            C4
blanchedalmond       FF           EB            CD
blue                 00           00            FF
blueviolet           8A           2B            E2
brown                A5           2A            2A
burlywood            DE           B8            87
cadetblue            5F           9E            A0
chartreuse           7F           FF            00
chocolate            D2           69            1E
coral                FF           7F            50
cornflowerblue       64           95            ED
cornsilk             FF           F8            DC
crimson              DC           14            3C
cyan                 00           FF            FF
darkblue             00           00            8B
darkcyan             00           8B            8B
darkgoldenrod        B8           86            0B
darkgray             A9           A9            A9
darkgreen            00           64            00
darkkhaki            BD           B7            6B
darkmagenta          8B           00            8B
darkolivegreen       55           6B            2F
darkorange           FF           8C            00
darkorchid           99           32            CC
darkred              8B           00            00
darksalmon           E9           96            7A
darkseagreen         8F           BC            8F
darkslateblue        48           3D            8B
darkslategray        2F           4F            4F
darkturquoise        00           CE            D1
darkviolet           94           00            D3
deeppink             FF           14            100
deepskyblue          00           BF            FF
dimgray              69           69            69
dodgerblue           1E           90            FF
firebrick            B2           22            22
floralwhite          FF           FA            F0
forestgreen          22           8B            22
fuchsia              FF           00            FF
gainsboro            DC           DC            DC
ghostwhite           F8           F8            FF
gold                 FF           D7            00
goldenrod            DA           A5            20
gray                 80           80            80
green                00           80            00
greenyellow          AD           FF            2F
honeydew             F0           FF            F0
hotpink              FF           69            B4
indianred            CD           5C            5C
indigo               4B           00            82
ivory                FF           FF            F0
khaki                F0           E6            8C
lavender             E6           E6            FA
lavenderblush        FF           F0            F5
lawngreen            7C           FC            00
lemonchiffon         FF           FA            CD
lightblue            AD           D8            E6
lightcoral           F0           80            80
lightcyan            E0           FF            FF
lightgoldenrodyellow FA           FA            D2
lightgreen           90           EE            90
lightgrey            D3           D3            D3
lightpink            FF           B6            C1
lightsalmon          FF           A0            7A
lightseagreen        20           B2            AA
lightskyblue         87           CE            FA
lightslategray       77           88            99
lightsteelblue       B0           C4            DE
lightyellow          FF           FF            E0
lime                 00           FF            00
limegreen            32           CD            32
linen                FA           F0            E6
magenta              FF           00            FF
maroon               80           00            00
mediumaquamarine     66           CD            AA
mediumblue           00           00            CD
mediumorchid         BA           55            D3
mediumpurple         100          70            DB
mediumseagreen       3C           B3            71
mediumslateblue      7B           68            EE
mediumspringgreen    00           FA            9A
mediumturquoise      48           D1            CC
mediumvioletred      C7           15            85
midnightblue         19           19            70
mintcream            F5           FF            FA
mistyrose            FF           E4            E1
moccasin             FF           E4            B5
navajowhite          FF           DE            AD
navy                 00           00            80
oldlace              FD           F5            E6
olive                80           80            00
olivedrab            6B           8E            23
orange               FF           A5            00
orangered            FF           45            00
orchid               DA           70            D6
palegoldenrod        EE           E8            AA
palegreen            98           FB            98
paleturquoise        AF           EE            EE
palevioletred        DB           70            100
papayawhip           FF           EF            D5
peachpuff            FF           DA            B9
peru                 CD           85            3F
pink                 FF           C0            CB
plum                 DD           A0            DD
powderblue           B0           E0            E6
purple               80           00            80
red                  FF           00            00
rosybrown            BC           8F            8F
royalblue            41           69            E1
saddlebrown          8B           45            13
salmon               FA           80            72
sandybrown           F4           A4            60
seagreen             2E           8B            57
seashell             FF           F5            EE
sienna               A0           52            2D
silver               C0           C0            C0
skyblue              87           CE            EB
slateblue            6A           5A            CD
slategray            70           80            90
snow                 FF           FA            FA
springgreen          00           FF            7F
steelblue            46           82            B4
tan                  D2           B4            8C
teal                 00           80            80
thistle              D8           BF            D8
tomato               FF           63            47
turquoise            40           E0            D0
violet               EE           82            EE
wheat                F5           DE            B3
whitesmoke           F5           F5            F5
yellow               FF           FF            00
yellowgreen          9A           CD            32
__END__

=head1 NAME

Bio::Graphics::Panel - Generate GD images of Bio::Seq objects

=head1 SYNOPSIS

  use Bio::Graphics;
  use Bio::DB::BioFetch;  # or some other Bio::SeqI generator

  # get a Bio::SeqI object somehow
  my $bf     = Bio::DB::BioFetch->new;
  my $cosmid = $bf->getSeq_by_id('CEF58D5');

  my @features = $seq->all_SeqFeatures;
  my @CDS      = grep {$_->primary_tag eq 'CDS'}  @features;
  my @gene     = grep {$_->primary_tag eq 'gene'} @features;
  my @tRNAs    = grep {$_->primary_tag eq 'tRNA'} @features;

  # let the drawing begin...
  my $panel = Bio::Graphics::Panel->new(
				      -segment => $cosmid,
				      -width  => 800
				     );

  $panel->add_track(arrow => $cosmid,
	  	   -bump => 0,
		   -double=>1,
		   -tick => 2);

  $panel->add_track(transcript  => \@gene,
		   -bgcolor    =>  'blue',
		   -fgcolor    =>  'black',
		   -key        => 'Genes',
		   -bump       =>  +1,
		   -height     =>  10,
		   -label      => 1,
		   -description=> 1
		 ) ;

  $panel->add_track(transcript2  => \@CDS,
		    -bgcolor    =>  'cyan',
		    -fgcolor    =>  'black',
		    -key        => 'CDS',
		    -bump       =>  +1,
		    -height     =>  10,
		    -label      => \&cds_label,
		    -description=> \&cds_description,
		 );

  $panel->add_track(generic    => \@tRNAs,
		    -bgcolor   =>  'red',
		    -fgcolor   =>  'black',
		    -key       => 'tRNAs',
		    -bump      =>  +1,
		    -height    =>  8,
		    -label      => 1,
		   );

  my $gd = $panel->gd;
  print $gd->can('png') ? $gd->png : $gd->gif;

  # these are callbacks used to generate nice labels and descriptions for
  # the features...
  sub cds_label {
    my $feature = shift;
    my @notes;
    foreach (qw(product gene)) {
      next unless $feature->has_tag($_);
      @notes = $feature->each_tag_value($_);
      last;
    }
    $notes[0];
  }

  sub cds_description {
    my $feature = shift;
    my @notes = $feature->each_tag_value('notes')
                if $feature->has_tag('notes');
    return unless @notes;
    substr($notes[0],30) = '...' if length $notes[0] > 30;
    $notes[0];
  }

=head1 DESCRIPTION

The Bio::Graphics::Panel class provides drawing and formatting
services for any object that implements the Bio::SeqFeatureI
interface, including Ace::Sequence::Feature and Das::Segment::Feature
objects.  It can be used to draw sequence annotations, physical
(contig) maps, or any other type of map in which a set of discrete
ranges need to be laid out on the number line.

The module supports a drawing style in which each type of feature
occupies a discrete "track" that spans the width of the display.  Each
track will have its own distinctive "glyph", a configurable graphical
representation of the feature.

The module also supports a more flexible style in which several
different feature types and their associated glyphs can occupy the
same track.  The choice of glyph is under run-time control.

Semantic zooming (for instance, changing the type of glyph depending
on the density of features) is supported by a callback system for
configuration variables.  The module has built-in support for Bio::Das
stylesheets, and stylesheet-driven configuration can be intermixed
with semantic zooming, if desired.

You can add a key to the generated image using either of two key
styles.  One style places the key captions at the top of each track.
The other style generates a graphical key at the bottom of the image.

Note that this modules depends on GD.

=head1 METHODS

This section describes the class and object methods for
Bio::Graphics::Panel.

Typically you will begin by creating a new Bio::Graphics::Panel
object, passing it the desired width of the image to generate and an
origin and length describing the coordinate range to display.  The
Bio::Graphics::Panel-E<gt>new() method has may configuration variables
that allow you to control the appearance of the image.

You will then call add_track() one or more times to add sets of
related features to the picture.  add_track() places a new horizontal
track on the image, and is likewise highly configurable.  When you
have added all the features you desire, you may call png() to convert
the image into a PNG-format image, or boxes() to return coordinate
information that can be used to create an imagemap.

=head2 CONSTRUCTORS

new() is the constructor for Bio::Graphics::Panel:

=over 4

=item $panel = Bio::Graphics::Panel-E<gt>new(@options)

The new() method creates a new panel object.  The options are
a set of tag/value pairs as follows:

  Option      Value                                  Default
  ------      -----                                  -------

  -offset     Base pair to place at extreme left     none
	      of image, in zero-based coordinates

  -length     Length of sequence segment, in bp      none

  -start      Start of range, in 1-based             none
              coordinates.

  -stop       Stop of range, in 1-based              none
	      coordinates.

  -segment    A Bio::SeqI or Das::Segment            none
              object, used to derive sequence
	      range if not otherwise specified.

  -width      Desired width of image, in pixels      600

  -spacing    Spacing between tracks, in pixels      5

  -pad_top    Additional whitespace between top      0
	      of image and contents, in pixels

  -pad_bottom Additional whitespace between top      0
	      of image and bottom, in pixels

  -pad_left   Additional whitespace between left     0
	      of image and contents, in pixels

  -pad_right  Additional whitespace between right    0
	      of image and bottom, in pixels

  -bgcolor    Background color for the panel as a    white
	      whole

  -key_color  Background color for the key printed   wheat
              at bottom of panel (if any)

  -key_spacing Spacing between key glyphs in the     10
               key printed at bottom of panel
               (if any)

  -key_font    Font to use in printed key            gdMediumBoldFont
	       captions.

  -key_style   Whether to print key at bottom of     none
	       panel ("bottom"), between each
	       track ("between"), or not at all
	       ("none").

  -all_callbacks Whether to invoke callbacks on      false
               the automatic "track" and "group"
               glyphs.

  -grid        Whether to draw a vertical grid in    false
               the background.  Pass a scalar true
               value to have a grid drawn at
               regular intervals (corresponding
               to the minor ticks of the arrow
	       glyph).  Pass an array reference
               to draw the grid at the specified
               positions.

  -gridcolor   Color of the grid                     lightcyan

Typically you will pass new() an object that implements the
Bio::RangeI interface, providing a length() method, from which the
panel will derive its scale.

  $panel = Bio::Graphics::Panel->new(-segment => $sequence,
				     -width   => 800);

new() will return undef in case of an error.

=back

=head2 OBJECT METHODS

=over 4

=item $track = $panel-E<gt>add_track($glyph,$features,@options)

The add_track() method adds a new track to the image. 

Tracks are horizontal bands which span the entire width of the panel.
Each track contains a number of graphical elements called "glyphs",
corresponding to a sequence feature. 

There are a large number of glyph types.  By default, each track will
be homogeneous on a single glyph type, but you can mix several glyph
types on the same track by providing a code reference to the -glyph
argument.  Other options passed to add_track() control the color and
size of the glyphs, whether they are allowed to overlap, and other
formatting attributes.  The height of a track is determined from its
contents and cannot be directly influenced.

The first two arguments are the glyph name and an array reference
containing the list of features to display.  The order of the
arguments is irrelevant, allowing either of these idioms:

  $panel->add_track(arrow => \@features);
  $panel->add_track(\@features => 'arrow');


The glyph name indicates how each feature is to be rendered.  A
variety of glyphs are available, and the number is growing. You may
omit the glyph name entirely by providing a B<-glyph> argument among
@options, as described below.

Currently, the following glyphs are available:

  Name        Description
  ----        -----------

  arrow	      An arrow; can be unidirectional or bidirectional.
	      It is also capable of displaying a scale with
	      major and minor tickmarks, and can be oriented
	      horizontally or vertically.

  cds         Draws CDS features, using the phase information to
              show the reading frame usage.  At high magnifications
              draws the protein translation.

  crossbox    A box with a big "X" inside it.

  diamond     A diamond, useful for point features like SNPs.

  dna         At high magnification draws the DNA sequence.  At
              low magnifications draws the GC content.

  dot         A circle, useful for point features like SNPs, stop
              codons, or promoter elements.

  ellipse     An oval.

  extending_arrow
              Similar to arrow, but a dotted line indicates when the
              feature extends beyond the end of the canvas.

  generic     A filled rectangle, nondirectional.

  graded_segments
              Similar to segments, but the intensity of the color
              is proportional to the score of the feature.  This
              is used for showing the intensity of blast hits or
              other alignment features.

  group	      A group of related features connected by a dashed line.
	      This is used internally by Panel.

  heterogeneous_segments
              Like segments, but you can use the source field of the feature
              to change the color of each segment.

  line        A simple line.

  pinsertion  A triangle designed to look like an insertion location
              (e.g. a transposon insertion).

  primers     Two inward pointing arrows connected by a line.
	      Used for STSs.

  rndrect     A round-cornered rectangle.

  segments    A set of filled rectangles connected by solid lines.
	      Used for interrupted features, such as gapped
	      alignments.

  ruler_arrow An arrow with major and minor tick marks and interval
              labels.

  toomany     Tries to show many features as a cloud.  Not very successful.

  track	      A group of related features not connected by a line.
	      This is used internally by Panel.

  transcript  Similar to segments, but the connecting line is
	      a "hat" shape, and the direction of transcription
	      is indicated by a small arrow.

  transcript2  Similar to transcript, but the direction of
              transcription is indicated by a terminal exon
              in the shape of an arrow.

  translation 1, 2 and 3-frame translations.  At low magnifications,
              can be configured to show start and stop codon locations.
              At high magnifications, shows the multi-frame protein
              translation.

  triangle    A triangle whose width and orientation can be altered.

If the glyph name is omitted from add_track(), the "generic" glyph
will be used by default.  To get more information about a glyph, run
perldoc on "Bio::Graphics::Glyph::glyphname", replacing "glyphname"
with the name of the glyph you are interested in.

The @options array is a list of name/value pairs that control the
attributes of the track.  Some options are interpretered directly by
the track.  Others are passed down to the individual glyphs (see
L<"GLYPH OPTIONS">).  The following options are track-specific:

  Option      Description                  Default
  ------      -----------                  -------

  -tkcolor    Track color                  white

  -glyph      Glyph class to use.         "generic"

  -stylesheet Bio::Das::Stylesheet to     none
              use to generate glyph
	      classes and options.

B<-tkcolor> controls the background color of the track as a whole.

B<-glyph> controls the glyph type.  If present, it supersedes the
glyph name given in the first or second argument to add_track().  The
value of B<-glyph> may be a constant string, a hash reference, or a
code reference.  In the case of a constant string, that string will be
used as the class name for all generated glyphs.  If a hash reference
is passed, then the feature's primary_tag() will be used as the key to
the hash, and the value, if any, used to generate the glyph type.  If
a code reference is passed, then this callback will be passed each
feature in turn as its single argument.  The callback is expected to
examine the feature and return a glyph name as its single result.

Example:

  $panel->add_track(\@exons,
		    -glyph => sub { my $feature = shift;
                                    $feature->source_tag eq 'curated'
                                          ? 'ellipse' : 'generic'; }
                    );

The B<-stylesheet> argument is used to pass a Bio::Das stylesheet
object to the panel.  This stylesheet will be called to determine both
the glyph and the glyph options.  If both a stylesheet and direct
options are provided, the latter take precedence.

If successful, add_track() returns an Bio::Graphics::Glyph object.
You can use this object to add additional features or to control the
appearance of the track with greater detail, or just ignore it.
Tracks are added in order from the top of the image to the bottom.  To
add tracks to the top of the image, use unshift_track().

B<Adding groups of features:> It is not uncommon to add a group of
features which are logically connected, such as the 5' and 3' ends of
EST reads.  To group features into sets that remain on the same
horizontal position and bump together, pass the sets as an anonymous
array.  For example:

  $panel->add_track(segments => [[$abc_5,$abc_3],
				 [$xxx_5,$xxx_3],
				 [$yyy_5,$yyy_3]]
		    );

Typical usage is:

 $panel->add_track( transcript    => \@genes,
 		    -fillcolor =>  'green',
 		    -fgcolor   =>  'black',
 		    -bump      =>  +1,
 		    -height    => 10,
 		    -label     => 1);

=item $track = unshift_track($glyph,$features,@options)

unshift_track() works like add_track(), except that the new track is
added to the top of the image rather than the bottom.

=item $gd = $panel-E<gt>gd

The gd() method lays out the image and returns a GD::Image object
containing it.  You may then call the GD::Image object's png() or
jpeg() methods to get the image data.

=item $png = $panel-E<gt>png

The png() method returns the image as a PNG-format drawing, without
the intermediate step of returning a GD::Image object.

=item $boxes = $panel-E<gt>boxes

=item @boxes = $panel-E<gt>boxes

The boxes() method returns the coordinates of each glyph, useful for
constructing an image map.  In a scalar context, boxes() returns an
array ref.  In an list context, the method returns the array directly.

Each member of the list is an anonymous array of the following format:

  [ $feature, $x1, $y1, $x2, $y2 ]

The first element is the feature object; either an
Ace::Sequence::Feature, a Das::Segment::Feature, or another Bioperl
Bio::SeqFeatureI object.  The coordinates are the topleft and
bottomright corners of the glyph, including any space allocated for
labels.

=back

=head1 GLYPH OPTIONS

Each glyph has its own specialized subset of options, but
some are shared by all glyphs:

  Option      Description                  Default
  ------      -----------                  -------

  -fgcolor    Foreground color		   black

  -bgcolor    Background color             turquoise

  -linewidth  Width of lines drawn by	   1
	      glyph

  -height     Height of glyph		   10

  -font       Glyph font		   gdSmallFont

  -fontcolor  Primary font color	   black

  -font2color Secondary font color	   turquoise

  -label      Whether to draw a label	   false

  -description  Whether to draw a          false
              description

  -bump	      Bump direction		   0

  -connector  Type of connector to         none
	      use to connect related
	      features.  Options are
	      "hat", "dashed" and
	      "none".

  -key        Description of track for     undef
	      use in key.

  -all_callbacks Whether to invoke         undef
              callbacks for autogenerated
              "track" and "group" glyphs

B<Specifying colors:> Colors can be expressed in either of two ways:
as symbolic names such as "cyan" and as HTML-style #RRGGBB triples.
The symbolic names are the 140 colors defined in the Netscape/Internet
Explorer color cube, and can be retrieved using the
Bio::Graphics::Panel-E<gt>color_names() method.

B<Foreground color:> The -fgcolor option controls the foreground
color, including the edges of boxes and the like.

B<Background color:> The -bgcolor option controls the background used
for filled boxes and other "solid" glyphs.  The foreground color
controls the color of lines and strings.  The -tkcolor argument
controls the background color of the entire track.

B<Track color:> The -tkcolor option used to specify the background of
the entire track.

B<Font color:> The -fontcolor option controls the color of primary
text, such as labels

B<Secondary Font color:> The -font2color option controls the color of
secondary text, such as descriptions.

B<Labels:> The -label argument controls whether or not the ID of the
feature should be printed next to the feature.  It is accepted by all
glyphs.  By default, the label is printed just above the glyph and
left aligned with it.  

-label can be a constant string or a code reference.  Values can be
any of:

  -label value     Description
  ------------     -----------

    0              Don't draw a label
    1              Calculate a label based on primary tag of sequence
    "a string"     Use "a string" as the label
    code ref       Invoke the code reference to compute the label

A known bug with this naming scheme is that you can't label a feature
with the string "1".  To work around this, use "1 " (note the terminal 
space).

B<Descriptions:> The -description argument controls whether or not a
brief description of the feature should be printed next to it.  By
default, the description is printed just below the glyph and
left-aligned with it.  A value of 0 will suppress the description.  A
value of 1 will call the source_tag() method of the feature.  A code
reference will be invoked to calculate the description on the fly.
Anything else will be treated as a string and used verbatim.

B<Connectors:> A glyph can contain subglyphs, recursively.  The top
level glyph is the track, which contains one or more groups, which
contain features, which contain subfeatures, and so forth.  By
default, the "group" glyph draws dotted lines between each of its
subglyphs, the "segment" glyph draws a solid line between each of its
subglyphs, and the "transcript" and "transcript2" glyphs draw
hat-shaped lines between their subglyphs.  All other glyphs do not
connect their components.  You can override this behavior by providing 
a -connector option, to explicitly set the type of connector.  Valid
options are "dashed", "solid", "hat" and "none".

B<Collision control:> The -bump argument controls what happens when
glyphs collide.  By default, they will simply overlap (value 0).  A
-bump value of +1 will cause overlapping glyphs to bump downwards
until there is room for them.  A -bump value of -1 will cause
overlapping glyphs to bump upwards.  Bump values of +2 and -2
implement a simpler bump algorithm in which each horizontal position
is occupied by one and only one feature.  The bump argument can also
be a code reference; see below.

B<Keys:> The -key argument declares that the track is to be shown in a
key appended to the bottom of the image.  The key contains a picture
of a glyph and a label describing what the glyph means.  The label is
specified in the argument to -key.

=head2 Options and Callbacks

Instead of providing a constant value to an option, you may subsitute
a code reference.  This code reference will be called every time the
panel needs to configure a glyph.  The callback will be called with
three arguments like this:

   sub callback {
      my ($feature,$option_name,$part_no,$total_parts,$glyph) = @_;
      # do something which results in $option_value being set
      return $option_value;
   }

The five arguments are C<$feature>, a reference to the IO::SeqFeatureI
object, C<$option_name>, the name of the option to configure,
C<$part_no>, an integer index indicating which subpart of the feature
is being drawn, C<$total_parts>, an integer indicating the total
number of subfeatures in the feature, and finally C<$glyph>, the Glyph
object itself.  The latter fields are useful in the common case of
treating the first or last subfeature differently, such as using a
different color for the terminal exon of a gene.

The callback should return a string indicating the desired value of
the option.  To tell the panel to use the default value for this
option, return the string "*default*".

When you install a callback for a feature that contains subparts, the
callback will be invoked first for the top-level feature, and then for
each of its subparts (recursively).  You should make sure to examine
the feature's type to determine whether the option is appropriate.

Some glyphs deliberately disable this recursive feature.  The "track",
"group", "transcript", "transcript2" and "segments" glyphs selectively
disable the -bump, -label and -description options.  This is to avoid,
for example, a label being attached to each exon in a transcript, or
the various segments of a gapped alignment bumping each other.  You
can override this behavior and force your callback to be invoked by
providing add_track() with a true B<-all_callbacks> argument.  In this 
case, you must be prepared to handle configuring options for the
"group" and "track" glyphs.

=head2 ACCESSORS

The following accessor methods provide access to various attributes of
the panel object.  Called with no arguments, they each return the
current value of the attribute.  Called with a single argument, they
set the attribute and return its previous value.

Note that in most cases you must change attributes prior to invoking
gd(), png() or boxes().  These three methods all invoke an internal
layout() method which places the tracks and the glyphs within them,
and then caches the result.

   Accessor Name      Description
   -------------      -----------

   width()	      Get/set width of panel
   spacing()	      Get/set spacing between tracks
   key_spacing()      Get/set spacing between keys
   length()	      Get/set length of segment (bp)
   pad_top()	      Get/set top padding
   pad_left()	      Get/set left padding
   pad_bottom()	      Get/set bottom padding
   pad_right()	      Get/set right padding

=head2 COLOR METHODS

The following methods are used internally, but may be useful for those
implementing new glyph types.

=over 4

=item @names = Bio::Graphics::Panel-E<gt>color_names

Return the symbolic names of the colors recognized by the panel
object.  In a scalar context, returns an array reference.

=item ($red,$green,$blue) = Bio::Graphics::Panel-E<gt>color_name_to_rgb($color)

Given a symbolic color name, returns the red, green, blue components
of the color.  In a scalar context, returns an array reference to the
rgb triplet.  Returns undef for an invalid color name.

=item @rgb = $panel-E<gt>rgb($index)

Given a GD color index (between 0 and 140), returns the RGB triplet
corresponding to this index.  This method is only useful within a
glyph's draw() routine, after the panel has allocated a GD::Image and
is populating it.

=item $index = $panel-E<gt>translate_color($color)

Given a color, returns the GD::Image index.  The color may be
symbolic, such as "turquoise", or a #RRGGBB triple, as in #F0E0A8.
This method is only useful within a glyph's draw() routine, after the
panel has allocated a GD::Image and is populating it.

=item $panel-E<gt>set_pen($width,$color)

Changes the width and color of the GD drawing pen to the values
indicated.  This is called automatically by the GlyphFactory fgcolor()
method.  It returns the GD value gdBrushed, which should be used for
drawing.

=back

=head1 BUGS

Please report them.

=head1 SEE ALSO

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

