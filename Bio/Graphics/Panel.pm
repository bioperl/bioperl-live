package Bio::Graphics::Panel;

use strict;
use Bio::Graphics::Glyph::Factory;
use Bio::Graphics::Feature;

# KEYLABELFONT must be treated as string until image_class is established
use constant KEYLABELFONT => 'gdMediumBoldFont';
use constant KEYSPACING   => 5; # extra space between key columns
use constant KEYPADTOP    => 5;  # extra padding before the key starts
use constant KEYCOLOR     => 'wheat';
use constant KEYSTYLE     => 'bottom';
use constant KEYALIGN     => 'left';
use constant GRIDCOLOR    => 'lightcyan';
use constant MISSING_TRACK_COLOR =>'gray';
use constant EXTRA_RIGHT_PADDING => 30;

use base qw(Bio::Root::Root);

my %COLORS;  # translation table for symbolic color names to RGB triple
my $IMAGEMAP = 'bgmap00001';
read_colors();

sub api_version { 1.654 }

# Create a new panel of a given width and height, and add lists of features
# one by one
sub new {
  my $class = shift;
  $class    = ref($class) || $class;
  my %options = @_;

  $class->read_colors() unless %COLORS;

  my $length = $options{-length} || 0;
  my $offset = $options{-offset}  || 0;
  my $spacing = $options{-spacing} || 5;
  my $bgcolor = $options{-bgcolor} || 'white';
  my $keyfont = $options{-key_font} || KEYLABELFONT;
  my $keycolor = $options{-key_color} || KEYCOLOR;
  my $keyspacing = $options{-key_spacing} || KEYSPACING;
  my $keystyle = $options{-key_style} || KEYSTYLE;
  my $keyalign = $options{-key_align} || KEYALIGN;
  my $allcallbacks = $options{-all_callbacks} || 0;
  my $gridcolor    = $options{-gridcolor} || GRIDCOLOR;
  my $grid         = $options{-grid}       || 0;
  my $extend_grid  = $options{-extend_grid}|| 0;
  my $flip         = $options{-flip}       || 0;
  my $empty_track_style   = $options{-empty_tracks} || 'key';
  my $autopad      = defined $options{-auto_pad} ? $options{-auto_pad} : 1;
  my $truecolor    = $options{-truecolor}  || 0;
  my $image_class  = ($options{-image_class} && $options{-image_class} =~ /SVG/)
                      ? 'GD::SVG'
		      : $options{-image_class} || 'GD';  # Allow users to specify GD::SVG using SVG
  my $linkrule     = $options{-link};
  my $titlerule    = $options{-title};
  my $targetrule   = $options{-target};
  my $background   = $options{-background};
  my $postgrid     = $options{-postgrid};
  $options{-stop}||= $options{-end};  # damn damn damn
  my $add_categories= $options{-add_category_labels};

  if (my $seg = $options{-segment}) {
    $offset = eval {$seg->start-1} || 0;
    $length = $seg->length;
  }

  $offset   ||= $options{-start}-1 if defined $options{-start};
  $length   ||= $options{-stop}-$options{-start}+1 
     if defined $options{-start} && defined $options{-stop};

  # bring in the image generator class, since we will need it soon anyway
  eval "require $image_class; 1" or $class->throw($@);

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
		extend_grid    => $extend_grid,
		bgcolor => $bgcolor,
		height => 0, # AUTO
		spacing => $spacing,
		key_font => $keyfont,
		key_color => $keycolor,
		key_spacing => $keyspacing,
		key_style => $keystyle,
		key_align => $keyalign,
		background => $background,
		postgrid   => $postgrid,
		autopad   => $autopad,
		all_callbacks => $allcallbacks,
		truecolor     => $truecolor,
		flip          => $flip,
		linkrule      => $linkrule,
		titlerule     => $titlerule,
		targetrule    => $targetrule,
		empty_track_style  => $empty_track_style,
		image_class  => $image_class,
		image_package => $image_class . '::Image',     # Accessors
		polygon_package => $image_class . '::Polygon',
		add_category_labels => $add_categories,
		key_boxes  => [],
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
sub extend_grid {
  my $self = shift;
  my $g = $self->{extend_grid};
  $self->{extend_grid} = shift if @_;
  $g;
}
sub flip {
  my $self = shift;
  my $g = $self->{flip};
  $self->{flip} = shift if @_;
  $g;
}

# values of empty_track_style are:
#    "suppress" -- suppress empty tracks entirely (default)
#    "key"      -- show just the key in "between" mode
#    "line"     -- draw a thin grey line
#    "dashed"   -- draw a dashed line
sub empty_track_style {
  my $self = shift;
  my $g = $self->{empty_track_style};
  $self->{empty_track_style} = shift if @_;
  $g;
}

sub key_style {
  my $self = shift;
  my $g = $self->{key_style};
  $self->{key_style} = shift if @_;
  $g;
}

sub auto_pad {
  my $self = shift;
  my $g = $self->{autopad};
  $self->{autopad} = shift if @_;
  $g;
}

# public routine for mapping from a base pair
# location to pixel coordinates
sub location2pixel {
  my $self   = shift;
  my $end    = $self->end + 1;
  my @coords = $self->{flip} ? map { $end-$_ } @_ : @_;
  $self->map_pt(@coords);
}

# numerous direct calls into array used here for performance considerations
sub map_pt {
  my $self   = shift;
  my $offset = $self->{offset};
  my $scale  = $self->{scale} || $self->scale;
  my $pl     = $self->{pad_left};
  my $pr     = $self->{width};
  my $flip   = $self->{flip};
  my $length = $self->{length};
  my @result;
  foreach (@_) {
    my $val = $flip 
      ? int (0.5 + $pr - ($length - ($_- 1)) * $scale)
      : int (0.5 + ($_-$offset-1) * $scale);
    $val = -1 if $val < 0;
    $val = $pr+1 if $val > $pr;
    push @result,$val;
  }
  @result;
}

sub map_no_trunc {
  my $self   = shift;
  my $offset = $self->{offset};
  my $scale  = $self->scale;
  my $pl     = $self->{pad_left};
  my $pr     = $pl + $self->{width}; # - $self->{pad_right};
  my $flip   = $self->{flip};
  my $length = $self->{length};
  my $end    = $offset+$length;
  my @result;
  foreach (@_) {
    my $val = $flip ? int (0.5 + $pl + ($end - ($_- 1)) * $scale) : int (0.5 + $pl + ($_-$offset-1) * $scale);
    push @result,$val;
  }
  @result;
}

sub scale {
  my $self = shift;
  # $self->{scale} ||= ($self->{width}-$self->pad_left-$self->pad_right)/($self->length);
  $self->{scale} ||= $self->width/($self->length);
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
  $self->pad_left + $self->width; # - $self->pad_right;
}
sub top {
  shift->pad_top;
}
sub bottom {
  my $self = shift;
  $self->height - $self->pad_bottom;
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

  local $^W = 0;  # uninitialized variable warnings under 5.00503

  my $panel_map =
    ref($map) eq 'CODE' ?  sub {
      my $feature = shift;
      return 'track' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'track' };
      return 'group' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'group' };
      return $map->($feature);
    }
   : ref($map) eq 'HASH' ? sub {
     my $feature = shift;
     return 'track' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'track' };
     return 'group' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'group' };
     return eval {$map->{$feature->primary_tag}} || 'generic';
   }
   : sub {
     my $feature = shift;
     return 'track' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'track' };
     return 'group' if eval { defined $feature->primary_tag && $feature->primary_tag  eq 'group' };
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
					    -start   => $self->offset+1,
					    -stop    => $self->offset+$self->length,
					    -type    => 'track'
					   );

  my $factory = Bio::Graphics::Glyph::Factory->new($self,@options);
  my $track   = $factory->make_glyph(-1,$feature);

  splice(@{$self->{tracks}},$position,0,$track);
  return $track;
}

sub _expand_padding {
  my $self   = shift;
  my $track  = shift;
  my $extra_padding = $self->extra_right_padding;

  my $keystyle          = $self->key_style;
  my $empty_track_style = $self->empty_track_style;

  return unless $keystyle eq 'left' or $keystyle eq 'right';
  return unless $self->auto_pad;

  $self->setup_fonts();
  my $width    = $self->{key_font}->width;

  my $key       = $self->track2key($track);
  return unless defined $key;

  my $has_parts = $track->parts;
  next if !$has_parts && $empty_track_style eq 'suppress';

  my $width_needed = $self->{key_font}->width * CORE::length($key)+3;
  if ($keystyle eq 'left') {
    my $width_i_have = $self->pad_left;
    $self->pad_left($width_needed)  if $width_needed > $width_i_have;
  } elsif ($keystyle eq 'right') {
    $width_needed += $extra_padding;
    my $width_i_have = $self->pad_right;
    $self->pad_right($width_needed) if $width_needed > $width_i_have;
  }
}

sub extra_right_padding { EXTRA_RIGHT_PADDING }

sub height {
  my $self = shift;
  $self->setup_fonts;

  for my $track (@{$self->{tracks}}) {
    $self->_expand_padding($track);
  }

  my $spacing           = $self->spacing;
  my $key_height        = $self->format_key;
  my $empty_track_style = $self->empty_track_style;
  my $key_style         = $self->key_style;
  my $bottom_key        = $key_style eq 'bottom';
  my $between_key       = $key_style eq 'between';
  my $side_key          = $key_style =~ /left|right/;
  my $draw_empty        = $empty_track_style =~ /^(line|dashed)$/;
  my $keyheight         = $self->{key_font}->height;
  my $height = 0;
  for my $track (@{$self->{tracks}}) {
    my $draw_between =  $between_key && $track->option('key');
    my $has_parts = $track->parts;
    next if !$has_parts && ($empty_track_style eq 'suppress'
		        or  $empty_track_style eq 'key' && $bottom_key);
    $height += $keyheight if $draw_between;
    $height += $self->spacing;
    my $layout_height = $track->layout_height;
    $height += ($side_key && $keyheight > $layout_height) ? $keyheight : $layout_height;
  }

  # get rid of spacing under last track
  $height -= $self->spacing unless $bottom_key;
  return $height + $key_height + $self->pad_top + $self->pad_bottom + 2;
}

sub setup_fonts {
  my $self = shift;
  return if ref $self->{key_font};

  my $image_class = $self->image_class;
  my $keyfont = $self->{key_font};
  my $font_obj = $image_class->$keyfont;
  $self->{key_font} = $font_obj;
}

sub gd {
  my $self        = shift;
  my $existing_gd = shift;

  local $^W = 0;  # can't track down the uninitialized variable warning

  return $self->{gd} if $self->{gd};

  $self->setup_fonts;

  unless ($existing_gd) {
    my $image_class = $self->image_class;
    eval "require $image_class; 1" or $self->throw($@);
  }

  my $height = $self->height;
  my $width  = $self->width + $self->pad_left + $self->pad_right;

  my $pkg = $self->image_package;
  my $gd  = $existing_gd || $pkg->new($width,$height,
				      ($self->{truecolor} && $pkg->can('isTrueColor') ? 1 : ())
				     );
  my %translation_table;
  for my $name ('white','black',keys %COLORS) {
    my $idx = $gd->colorAllocate(@{$COLORS{$name}});
    $translation_table{$name} = $idx;
  }

  $self->{translations} = \%translation_table;
  $self->{gd}           = $gd;
  if ($self->bgcolor) {
    $gd->fill(0,0,$self->bgcolor);
  } elsif (eval {$gd->isTrueColor}) {
    $gd->fill(0,0,$translation_table{'white'});
  }

  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $offset = $pt;
  my $keyheight   = $self->{key_font}->height;
  my $bottom_key  = $self->{key_style} eq 'bottom';
  my $between_key = $self->{key_style} eq 'between';
  my $left_key    = $self->{key_style} eq 'left';
  my $right_key   = $self->{key_style} eq 'right';
  my $empty_track_style = $self->empty_track_style;
  my $spacing = $self->spacing;

  # we draw in two steps, once for background of tracks, and once for
  # the contents.  This allows the grid to sit on top of the track background.
  for my $track (@{$self->{tracks}}) {
    my $draw_between = $between_key && $track->option('key');
    next if !$track->parts && ($empty_track_style eq 'suppress'
			   or  $empty_track_style eq 'key' && $bottom_key);
    $gd->filledRectangle($pl,
			 $offset,
			 $width-$self->pad_right,
			 $offset+$track->layout_height
			 + ($between_key ? $self->{key_font}->height : 0),
			 $track->tkcolor)
      if defined $track->tkcolor;
    $offset += $keyheight if $draw_between;
    $offset += $track->layout_height + $spacing;
  }

  $self->draw_background($gd,$self->{background})  if $self->{background};
  $self->draw_grid($gd)                            if $self->{grid};
  $self->draw_background($gd,$self->{postgrid})    if $self->{postgrid};

  $offset = $pt;
  for my $track (@{$self->{tracks}}) {
    my $draw_between = $between_key && $track->option('key');
    my $has_parts = $track->parts;
    my $side_key_height = 0;

    next if !$has_parts && ($empty_track_style eq 'suppress'
			or  $empty_track_style eq 'key' && $bottom_key);

    if ($draw_between) {
      $offset += $self->draw_between_key($gd,$track,$offset);
    }


    $self->draw_empty($gd,$offset,$empty_track_style)
      if !$has_parts && $empty_track_style=~/^(line|dashed)$/;

    $track->draw($gd,$pl,$offset,0,1);

    if ($self->{key_style} =~ /^(left|right)$/) {
      $side_key_height = $self->draw_side_key($gd,$track,$offset,$self->{key_style});
    }

    $self->track_position($track,$offset);
    my $layout_height = $track->layout_height;
    $offset += ($side_key_height > $layout_height ? $side_key_height : $layout_height)+$spacing;
  }


  $self->draw_bottom_key($gd,$pl,$offset) if $self->{key_style} eq 'bottom';
  return $self->{gd} = $gd;
}


# Package accessors
# GD (and GD::SVG)'s new() resides in GD::Image
sub image_class     { return shift->{image_class}; }
sub image_package   { return shift->{image_package}; }
sub polygon_package { return shift->{polygon_package}; }

sub boxes {
  my $self = shift;
  my @boxes;
  my $offset = 0;

  $self->setup_fonts;

  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $between_key       = $self->{key_style} eq 'between';
  my $bottom_key        = $self->{key_style} eq 'bottom';
  my $empty_track_style = $self->empty_track_style;
  my $keyheight         = $self->{key_font}->height;
  my $spacing = $self->spacing;

  for my $track (@{$self->{tracks}}) {
    my $draw_between =  $between_key && $track->option('key');
    next if !$track->parts && ($empty_track_style eq 'suppress'
			    or  $empty_track_style eq 'key' && $bottom_key);
    $offset += $keyheight if $draw_between;
    my $boxes = $track->boxes($pl,$offset+$pt);
    $self->track_position($track,$offset);
    push @boxes,@$boxes;
    $offset += $track->layout_height + $self->spacing;
  }
  return wantarray ? @boxes : \@boxes;
}

sub track_position {
  my $self  = shift;
  my $track = shift;
  my $d = $self->{_track_position}{$track};
  $self->{_track_position}{$track} = shift if @_;
  $d;
}

# draw the keys -- between
sub draw_between_key {
  my $self   = shift;
  my ($gd,$track,$offset) = @_;
  my $key = $self->track2key($track) or return 0;
  my $x =   $self->{key_align} eq 'center' ? $self->width - (CORE::length($key) * $self->{key_font}->width)/2
          : $self->{key_align} eq 'right'  ? $self->width - CORE::length($key)
          : $self->pad_left;

  # Key color hard-coded. Should be configurable for the control freaks.
  my $color = $self->translate_color('black');
  $gd->string($self->{key_font},$x,$offset,$key,$color);
  $self->add_key_box($track,$key,$x,$offset);
  return $self->{key_font}->height;
}

# draw the keys -- left or right side
sub draw_side_key {
  my $self   = shift;
  my ($gd,$track,$offset,$side) = @_;
  my $key = $self->track2key($track) or return;
  my $pos = $side eq 'left' ? $self->pad_left - $self->{key_font}->width * CORE::length($key)-3
                            : $self->pad_left + $self->width + EXTRA_RIGHT_PADDING;
  my $color = $self->translate_color('black');
  $gd->filledRectangle($pos,$offset,
		 $pos+$self->{key_font}->width*CORE::length($key),$offset,#-$self->{key_font}->height)/2,
		 $self->bgcolor);
  $gd->string($self->{key_font},$pos,$offset,$key,$color);
  $self->add_key_box($track,$key,$pos,$offset);
  return $self->{key_font}->height;
}

# draw the keys -- bottom
sub draw_bottom_key {
  my $self = shift;
  my ($gd,$left,$top) = @_;
  my $key_glyphs = $self->{key_glyphs} or return;

  my $color = $self->translate_color($self->{key_color});
  $gd->filledRectangle($left,$top,$self->width - $self->pad_right,$self->height-$self->pad_bottom,$color);
  my $text_color = $self->translate_color('black');
  $gd->string($self->{key_font},$left,KEYPADTOP+$top,"KEY:",$text_color);
  $top += $self->{key_font}->height + KEYPADTOP;
  $_->draw($gd,$left,$top) foreach @$key_glyphs;
}

# Format the key section, and return its height
sub format_key {
  my $self = shift;
  return 0 unless $self->key_style eq 'bottom';

  return $self->{key_height} if defined $self->{key_height};

  my $suppress = $self->{empty_track_style} eq 'suppress';
  my $between  = $self->{key_style}         eq 'between';

  if ($between) {
    my @key_tracks = $suppress
      ? grep {$_->option('key') && $_->parts} @{$self->{tracks}}
      : grep {$_->option('key')} @{$self->{tracks}};
    return $self->{key_height} = @key_tracks * $self->{key_font}->height;
  }

  elsif ($self->{key_style} eq 'bottom') {

    my ($height,$width) = (0,0);
    my %tracks;
    my @glyphs;
    local $self->{flip} = 0;  # don't want to worry about flipped keys!

    # determine how many glyphs become part of the key
    # and their max size
    for my $track (@{$self->{tracks}}) {

      next unless $track->option('key');
      next if $suppress && !$track->parts;

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

sub add_key_box {
  my $self = shift;
  my ($track,$label,$x,$y) = @_;
  my $value = [$label,$x,$y,$x+$self->{key_font}->width*CORE::length($label),$y+$self->{key_font}->height,$track];
  push @{$self->{key_boxes}},$value;
}

sub key_boxes {
  my $ref  = shift->{key_boxes};
  return wantarray ? @$ref : $ref;
}

sub add_category_labels {
  my $self = shift;
  my $d    = $self->{add_category_labels};
  $self->{add_category_labels} = shift if @_;
  $d;
}

sub track2key {
  my $self = shift;
  my $track = shift;
  return $track->make_key_name();
}

sub draw_empty {
  my $self  = shift;
  my ($gd,$offset,$style) = @_;
  $offset  += $self->spacing/2;
  my $left  = $self->pad_left;
  my $right = $self->width-$self->pad_right;
  my $color = $self->translate_color(MISSING_TRACK_COLOR);
  my $ic    = $self->image_class;
  if ($style eq 'dashed') {
    $gd->setStyle($color,$color,$ic->gdTransparent(),$ic->gdTransparent());
    $gd->line($left,$offset,$right,$offset,$ic->gdStyled());
  } else {
    $gd->line($left,$offset,$right,$offset,$color);
  }
  $offset;
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
    my $first_tick = $minor * int($self->start/$minor);
    for (my $i = $first_tick-1; $i <= $self->end+1; $i += $minor) {
      push @positions,$i;
    }
  }
  my $pl = $self->pad_left;
  my $pt = $self->extend_grid ? 0 : $self->pad_top;
  my $pr = $self->right;
  my $pb = $self->extend_grid ? $self->height : $self->height - $self->pad_bottom;
  my $offset = $self->{offset}+$self->{length}+1;
  for my $tick (@positions) {
    my ($pos) = $self->map_pt($self->{flip} ? $offset - $tick
                                            : $tick);

    $gd->line($pl+$pos,$pt,$pl+$pos,$pb,$gridcolor);
  }
}

# draw an image (or invoke a drawing routine)
sub draw_background {
  my $self = shift;
  my ($gd,$image_or_routine) = @_;
  if (ref $image_or_routine eq 'CODE') {
    return $image_or_routine->($gd,$self);
  }
  if (-f $image_or_routine) { # a file to draw
    my $method = $image_or_routine =~ /\.png$/i   ? 'newFromPng'
               : $image_or_routine =~ /\.jpe?g$/i ? 'newFromJpeg'
               : $image_or_routine =~ /\.gd$/i    ? 'newFromGd'
               : $image_or_routine =~ /\.gif$/i   ? 'newFromGif'
               : $image_or_routine =~ /\.xbm$/i   ? 'newFromXbm'
	       : '';
    return unless $method;
    my $image = eval {$self->image_package->$method($image_or_routine)};
    unless ($image) {
      warn $@;
      return;
    }
    my ($src_width,$src_height) = $image->getBounds;
    my ($dst_width,$dst_height) = $gd->getBounds;
    # tile the thing on
    for (my $x = 0; $x < $dst_width; $x += $src_width) {
      for (my $y = 0; $y < $dst_height; $y += $src_height) {
	$gd->copy($image,$x,$y,0,0,$src_width,$src_height);
      }
    }
  }
}

# calculate major and minor ticks, given a start position
sub ticks {
  my $self = shift;
  my ($length,$minwidth) = @_;

  my $img = $self->image_class;
  $length   = $self->{length}             unless defined $length;
  $minwidth = $img->gdSmallFont->width*7  unless defined $minwidth;

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

  # to make sure a major tick shows up somewhere in the first half
  #
  # $interval *= .5 if ($interval > 0.5*$length);

  return ($interval,$interval/10);
}

# reverse of translate(); given index, return rgb triplet
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
    return $self->colorClosest($gd,@colors);
  }
  elsif ($colors[0] =~ /^\#([0-9A-F]{2})([0-9A-F]{2})([0-9A-F]{2})$/i) {
    my $gd = $self->gd or return 1;
    my ($r,$g,$b) = (hex($1),hex($2),hex($3));
    return $self->colorClosest($gd,$r,$g,$b);
  }
  else {
    my $color = $colors[0];
    my $table = $self->{translations} or return 1;
    return defined $table->{$color} ? $table->{$color} : 1;
  }
}

# workaround for bad GD
sub colorClosest {
  my ($self,$gd,@c) = @_;
  return $self->{closestcache}{"@c"} if exists $self->{closestcache}{"@c"};
  return $self->{closestcache}{"@c"} = $gd->colorResolve(@c) if $GD::VERSION < 2.04;

  my $index = $gd->colorResolve(@c);
  return $self->{closestcache}{"@c"} = $index if $index >= 0;

  my $value;
  for (keys %COLORS) {
    my ($r,$g,$b) = @{$COLORS{$_}};
    my $dist = ($r-$c[0])**2 + ($g-$c[1])**2 + ($b-$c[2])**2;
    ($value,$index) = ($dist,$_) if !defined($value) || $dist < $value;
  }
  return $self->{closestcache}{"@c"} = $self->{translations}{$index};
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
  my $gd = $self->{gd};
  my $pkg = $self->image_package;
  my $pen = $self->{pens}{$linewidth} = $pkg->new($linewidth,$linewidth);
  my @rgb = $self->rgb($color);
  my $bg = $pen->colorAllocate(255,255,255);
  my $fg = $pen->colorAllocate(@rgb);
  $pen->fill(0,0,$fg);
  $gd->setBrush($pen);
  return $self->image_class->gdBrushed();
}

sub png {
  my $gd = shift->gd;
  $gd->png;
}

sub svg {
  my $gd = shift->gd;
  $gd->svg;
}


# WARNING: THIS STUFF IS COPIED FROM Bio::Graphics::Browser.pm AND
# Bio::Graphics::FeatureFile AND MUST BE REFACTORED
# write a png image to disk and generate an image map in a convenient
# CGIish way.
sub image_and_map {
  my $self        = shift;
  my %args        = @_;
  my $link_rule   = $args{-link}    || $self->{linkrule};
  my $title_rule  = $args{-title}   || $self->{titlerule};
  my $target_rule = $args{-target}  || $self->{targetrule};
  my $tmpurl      = $args{-url}     || '/tmp';
  my $docroot     = $args{-root}    || $ENV{DOCUMENT_ROOT} || '';
  my $mapname     = $args{-mapname} || $IMAGEMAP++;
  $docroot       .= '/' if $docroot && $docroot !~ m!/$!;

  # get rid of any netstat part please
  (my $tmpurlbase = $tmpurl) =~ s!^\w+://[^/]+!!;

  my $tmpdir    = "${docroot}${tmpurlbase}";

  my $url       = $self->create_web_image($tmpurl,$tmpdir);
  my $map       = $self->create_web_map($mapname,$link_rule,$title_rule,$target_rule);
  return ($url,$map,$mapname);
}

sub create_web_image {
  my $self             = shift;
  my ($tmpurl,$tmpdir) = @_;

  # create directory if it isn't there already
  # we need to untaint tmpdir before calling mkpath()
  return unless $tmpdir =~ /^(.+)$/;
  my $path = $1;
  unless (-d $path) {
    require File::Path unless defined &File::Path::mkpath;
    File::Path::mkpath($path,0,0777) or $self->throw("Couldn't create temporary image directory $path: $!");
  }

  unless (defined &Digest::MD5::md5_hex) {
    eval "require Digest::MD5; 1"
      or $self->throw("Sorry, but the image_and_map() method requires the Digest::MD5 module.");
  }
  my $data      = $self->png;
  my $signature = Digest::MD5::md5_hex($data);
  my $extension = 'png';

  # untaint signature for use in open
  $signature =~ /^([0-9A-Fa-f]+)$/g or return;
  $signature = $1;

  my $url         = sprintf("%s/%s.%s",$tmpurl,$signature,$extension);
  my $imagefile   = sprintf("%s/%s.%s",$tmpdir,$signature,$extension);

  open (my $F,">", $imagefile) || $self->throw("Can't open image file $imagefile for writing: $!\n");
  binmode($F);
  print $F $data;

  return $url;
}

sub create_web_map {
  my $self     = shift;
  my ($name,$linkrule,$titlerule,$targetrule) = @_;
  $name ||= 'map';
  my $boxes    = $self->boxes;
  my (%track2link,%track2title,%track2target);

  my $map = qq(<map name="$name" id="$name">\n);
  foreach (@$boxes){
    my ($feature,$left,$top,$right,$bottom,$track) = @$_;
    next unless $feature->can('primary_tag');

    my $lr  = $track2link{$track} ||= (defined $track->option('link') ? $track->option('link') : $linkrule);
    next unless   $lr;

    my $tr   = exists $track2title{$track} 
      ? $track2title{$track}
      : $track2title{$track} ||= (defined $track->option('title')  ? $track->option('title')  : $titlerule);
    my $tgr  = exists $track2target{$track} 
      ? $track2target{$track}
      : $track2target{$track} ||= (defined $track->option('target')? $track->option('target')  : $targetrule);

    my $href   = $self->make_link($lr,$feature);
    my $alt    = $self->make_link($tr,$feature);
    my $target = $self->make_link($tgr,$feature);
    $alt       = $self->make_title($feature) unless defined $alt;

    my $a      = $alt    ? qq(title="$alt" alt="$alt") : '';
    my $t      = $target ? qq(target="$target")        : '';
    $map .= qq(<area shape="rect" coords="$left,$top,$right,$bottom" href="$href" $a $t/>\n);
  }
  $map .= "</map>\n";
  $map;
}

sub make_link {
  my $self = shift;
  my ($linkrule,$feature) = @_;
  eval "require Bio::Graphics::FeatureFile;1"
    unless Bio::Graphics::FeatureFile->can('link_pattern');
  return Bio::Graphics::FeatureFile->link_pattern($linkrule,$feature,$self);
}

sub make_title {
  my $self = shift;
  my $feature = shift;
  eval "require Bio::Graphics::FeatureFile;1"
    unless Bio::Graphics::FeatureFile->can('make_title');
  return Bio::Graphics::FeatureFile->make_title($feature);
}

sub read_colors {
  my $class = shift;
  lock %COLORS;
  local ($/) = "\n";
  while (<DATA>) {
    chomp;
    last if /^__END__/;
    my ($name,$r,$g,$b) = split /\s+/;
    @{$COLORS{$name}} = (hex $r,hex $g, hex $b);
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

sub finished {
    my $self = shift;
    for my $track (@{$self->{tracks} || []}) {
	$track->finished();
    }
    delete $self->{tracks};
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
grey                 80           80            80
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
gradient1	00 ff 00
gradient2	0a ff 00
gradient3	14 ff 00
gradient4	1e ff 00
gradient5	28 ff 00
gradient6	32 ff 00
gradient7	3d ff 00
gradient8	47 ff 00
gradient9	51 ff 00
gradient10	5b ff 00
gradient11	65 ff 00
gradient12	70 ff 00
gradient13	7a ff 00
gradient14	84 ff 00
gradient15	8e ff 00
gradient16	99 ff 00
gradient17	a3 ff 00
gradient18	ad ff 00
gradient19	b7 ff 00
gradient20	c1 ff 00
gradient21	cc ff 00
gradient22	d6 ff 00
gradient23	e0 ff 00
gradient24	ea ff 00
gradient25	f4 ff 00
gradient26	ff ff 00
gradient27	ff f4 00
gradient28	ff ea 00
gradient29	ff e0 00
gradient30	ff d6 00
gradient31	ff cc 00
gradient32	ff c1 00
gradient33	ff b7 00
gradient34	ff ad 00
gradient35	ff a3 00
gradient36	ff 99 00
gradient37	ff 8e 00
gradient38	ff 84 00
gradient39	ff 7a 00
gradient40	ff 70 00
gradient41	ff 65 00
gradient42	ff 5b 00
gradient43	ff 51 00
gradient44	ff 47 00
gradient45	ff 3d 00
gradient46	ff 32 00
gradient47	ff 28 00
gradient48	ff 1e 00
gradient49	ff 14 00
gradient50	ff 0a 00
__END__

=head1 NAME

Bio::Graphics::Panel - Generate GD images of Bio::Seq objects

=head1 SYNOPSIS

 # This script parses a GenBank or EMBL file named on the command
 # line and produces a PNG rendering of it.  Call it like this:
 # render.pl my_file.embl | display -

 use strict;
 use Bio::Graphics;
 use Bio::SeqIO;

 my $file = shift                       or die "provide a sequence file as the argument";
 my $io = Bio::SeqIO->new(-file=>$file) or die "could not create Bio::SeqIO";
 my $seq = $io->next_seq                or die "could not find a sequence in the file";

 my @features = $seq->all_SeqFeatures;

 # sort features by their primary tags
 my %sorted_features;
 for my $f (@features) {
   my $tag = $f->primary_tag;
   push @{$sorted_features{$tag}},$f;
 }

 my $panel = Bio::Graphics::Panel->new(
                                      -length    => $seq->length,
 				      -key_style => 'between',
 				      -width     => 800,
 				      -pad_left  => 10,
 				      -pad_right => 10,
 				      );
 $panel->add_track( arrow => Bio::SeqFeature::Generic->new(-start=>1,
                                                           -end=>$seq->length),
 		  -bump => 0,
 		  -double=>1,
 		  -tick => 2);
 $panel->add_track(generic => Bio::SeqFeature::Generic->new(-start=>1,
							  -end=>$seq->length),
 		  -glyph  => 'generic',
 		  -bgcolor => 'blue',
 		  -label  => 1,
 		 );

 # general case
 my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
 my $idx    = 0;
 for my $tag (sort keys %sorted_features) {
   my $features = $sorted_features{$tag};
   $panel->add_track($features,
 		    -glyph    =>  'generic',
 		    -bgcolor  =>  $colors[$idx++ % @colors],
 		    -fgcolor  => 'black',
 		    -font2color => 'red',
 		    -key      => "${tag}s",
 		    -bump     => +1,
 		    -height   => 8,
 		    -label    => 1,
 		    -description => 1,
 		   );
 }

 print $panel->png;
 $panel->finished;

 exit 0;

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

Note that this module depends on GD. The optional SVG output depends
on GD::SVG and SVG.

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

  -end        Same as -stop.

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
	       track ("between"), to the left of
               each track ("left"), to the right
               of each track ("right") or
               not at all ("none").

  -add_category_labels                               false
               Whether to add the "category" to
               the track key. The category is
               an optional argument that can
               be attached to each track. If
               a category is present, and this
               option is true, then the category
               will be added to the track label
               in parentheses. For example, if
               -key is "Protein matches" and
               -category is "vertebrate", then
               the track will be labeled
               "Protein matches (vertebrate)".

  -auto_pad    If "left" or "right" keys are in use  true
               then setting auto_pad to a true value
               will allow the panel to adjust its
               width in order to accomodate the
               length of the longest key.

  -empty_tracks What to do when a track is empty.    suppress
              Options are to suppress the track
              completely ("suppress"), to show just
              the key in "between" mode ("key"),
              to draw a thin grey line ("line"),
              or to draw a dashed line ("dashed").

  -flip       flip the drawing coordinates left     false
              to right, so that lower coordinates
              are to the right.  This can be
              useful for drawing (-) strand
              features.

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

  -extend_grid If true, extend the grid into the pad false
               top and pad_bottom regions

  -background  An image or callback to use for the   none
               background of the image. Will be
               invoked I<before> drawing the grid.

  -postgrid    An image or callback to use for the   none
               background of the image.  Will be 
               invoked I<after> drawing the grid.

  -truecolor   Create a truecolor (24-bit) image.    false
               Useful when working with the
               "image" glyph.

  -image_class To create output in scalable vector
               graphics (SVG), optionally pass the image
               class parameter 'GD::SVG'. Defaults to
               using vanilla GD. See the corresponding
               image_class() method below for details.

  -link, -title, -target
               These options are used when creating imagemaps
               for display on the web.  See L</"Creating Imagemaps">.


Typically you will pass new() an object that implements the
Bio::RangeI interface, providing a length() method, from which the
panel will derive its scale.

  $panel = Bio::Graphics::Panel->new(-segment => $sequence,
				     -width   => 800);

new() will return undef in case of an error.

Note that if you use the "left" or "right" key styles, you are
responsible for allocating sufficient -pad_left or -pad_right room for
the labels to appear.  The necessary width is the number of characters
in the longest key times the font width (gdMediumBoldFont by default)
plus 3 pixels of internal padding.  The simplest way to calculate this
is to iterate over the possible track labels, find the largest one,
and then to compute its width using the formula:

  $width = gdMediumBoldFont->width * length($longest_key) +3;

In order to obtain scalable vector graphics (SVG) output, you should
pass new() the -image_class=E<gt>'GD::SVG' parameter. This will cause
Bio::Graphics::Panel to load the optional GD::SVG module. See the gd()
and svg() methods below for additional information.

You can tile an image onto the panel either before or after it draws
the grid. Simply provide the filename of the image in the -background
or -postgrid options. The image file must be of type PNG, JPEG, XBM or
GIF and have a filename ending in .png, .jpg, .jpeg, .xbm or .gif.

You can also pass a code ref for the -background or -postgrid option,
in which case the subroutine will be invoked at the appropriate time
with the GD::Image object and the Panel object as its two arguments.
You can then use the panel methods to map base pair coordinates into
pixel coordinates and do some custom drawing.  For example, this code
fragment will draw a gray rectangle between bases 500 and 600 to
indicate a "gap" in the sequence:

  my $panel = Bio::Graphics::Panel->new(-segment=>$segment,
                                        -grid=>1,
                                        -width=>600,
                                        -postgrid=> \&draw_gap);
  sub gap_it {
     my $gd    = shift;
     my $panel = shift;
     my ($gap_start,$gap_end) = $panel->location2pixel(500,600);
     my $top                  = $panel->top;
     my $bottom               = $panel->bottom;
     my $gray                 = $panel->translate_color('gray');
     $gd->filledRectangle($gap_start,$top,$gap_end,$bottom,$gray);
}


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

  anchored_arrow
              a span with vertical bases |---------|.  If one or
              the other end of the feature is off-screen, the base
              will be replaced by an arrow.

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

  image	      A pixmap image that will be layered on top of the graphic.

  heterogeneous_segments
              Like segments, but you can use the source field of the feature
              to change the color of each segment.

  line        A simple line.

  pinsertion  A triangle designed to look like an insertion location
              (e.g. a transposon insertion).

  processed_transcript  multi-purpose representation of a spliced mRNA, including
			positions of UTRs

  primers     Two inward pointing arrows connected by a line.
	      Used for STSs.

  redgreen_box A box that changes from green->yellow->red as the score
              of the feature increases from 0.0 to 1.0.  Useful for
              representing microarray results.

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

  xyplot      Histograms and other graphs plotted against the genome.

  whiskerplot Box and whisker plot for statistical data

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

=item $track = $panel-E<gt>insert_track($position,$glyph,$features,@options)

This works like add_track(), but the track is inserted into the
indicated position.  The track will be inserted B<before> the
indicated position; thus specify a track of 0 to insert the new track
at the beginning.

=item $gd = $panel-E<gt>gd([$gd])

The gd() method lays out the image and returns a GD::Image object
containing it.  You may then call the GD::Image object's png() or
jpeg() methods to get the image data.

Optionally, you may pass gd() a preexisting GD::Image object that you
wish to draw on top of.  If you do so, you should call the width() and
height() methods first to ensure that the image has sufficient
dimensions.

If you passed new() the -image_class=E<gt>'GD::SVG' parameter, the gd() method
returns a GD::SVG::Image object. This object overrides GD::Image
methods in order to generate SVG output. It behaves exactly as
described for GD::Image objects with one exception: it implements and
svg() method instead of the png() or jpeg() methods. Currently there
is no direct access to underlying SVG calls but this is subject to
change in the future.

=item $png = $panel-E<gt>png

The png() method returns the image as a PNG-format drawing, without
the intermediate step of returning a GD::Image object.

=item $svg = $panel-E<gt>svg

The svg() method returns the image in an XML-ified SVG format.

=item $panel-E<gt>finished

Bio::Graphics creates memory cycles.  When you are finished with the
panel, you should call its finished() method.  Otherwise you will have
memory leaks.  This is only an issue if you're going to create several
panels in a single program.

=item $image_class = $panel-E<gt>image_class

The image_class() method returns the current drawing package being
used, currently one of GD or GD::SVG.  This is primarily used
internally to ensure that calls to GD's exported methods are called in
an object-oriented manner to avoid compile time undefined string
errors.  This is usually not needed for external use.

=item $image_package = $panel-E<gt>image_package

This accessor method, like image_class() above is provided as a
convenience.  It returns the current image package in use, currently
one of GD::Image or GD::SVG::Image.  This is not normally used
externally.

=item $polygon_package = $panel-E<gt>polygon_package

This accessor method, like image_package() above is provided as a
convenience.  It returns the current polygon package in use, currently
one of GD::Polygon or GD::SVG::Polygon.  This is not normally used
externally except in the design of glyphs.

=item $boxes = $panel-E<gt>boxes

=item @boxes = $panel-E<gt>boxes

The boxes() method returns a list of arrayrefs containing the
coordinates of each glyph.  The method is useful for constructing an
image map.  In a scalar context, boxes() returns an arrayref.  In an
list context, the method returns the list directly.

Each member of the list is an arrayref of the following format:

  [ $feature, $x1, $y1, $x2, $y2, $track ]

The first element is the feature object; either an
Ace::Sequence::Feature, a Das::Segment::Feature, or another Bioperl
Bio::SeqFeatureI object.  The coordinates are the topleft and
bottomright corners of the glyph, including any space allocated for
labels. The track is the Bio::Graphics::Glyph object corresponding to
the track that the feature is rendered inside.

=item $boxes = $panel-E<gt>key_boxes

=item @boxes = $panel-E<gt>key_boxes

Returns the positions of the track keys as an arrayref or a list,
depending on context. Each value in the list is an arrayref of format:

 [ $key_text, $x1, $y1, $x2, $y2, $track ]

=item $position = $panel-E<gt>track_position($track)

After calling gd() or boxes(), you can learn the resulting Y
coordinate of a track by calling track_position() with the value
returned by add_track() or unshift_track().  This will return undef if
called before gd() or boxes() or with an invalid track.

=item @pixel_coords = $panel-E<gt>location2pixel(@feature_coords)

Public routine to map feature coordinates (in base pairs) into pixel
coordinates relative to the left-hand edge of the picture. If you
define a -background callback, the callback may wish to invoke this
routine in order to translate base coordinates into pixel coordinates.

=item $left = $panel-E<gt>left

=item $right = $panel-E<gt>right

=item $top   = $panel-E<gt>top

=item $bottom = $panel-E<gt>bottom

Return the pixel coordinates of the I<drawing area> of the panel, that
is, exclusive of the padding.

=back

=head1 GLYPH OPTIONS

Each glyph has its own specialized subset of options, but
some are shared by all glyphs:

  Option      Description                  Default
  ------      -----------                  -------

  -key        Description of track for     undef
	      display in the track label.

  -category   The category of the track    undef
	      for display in the
              track label.

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

  -sort_order Specify layout sort order    "default"

  -bump_limit Maximum number of levels     undef (unlimited)
              to bump

  -hbumppad   Additional horizontal        0
              padding between bumped
              features

  -strand_arrow Whether to indicate        undef (false)
                 strandedness

  -stranded    Synonym for -strand_arrow   undef (false)

  -part_labels Whether to label individual undef (false)
               subparts.

  -part_label_merge Whether to merge       undef (false)
              adjacent subparts when
              labeling.

  -connector  Type of connector to         none
	      use to connect related
	      features.  Options are
	      "solid," "hat", "dashed", 
              "quill" and "none".

  -all_callbacks Whether to invoke         undef
              callbacks for autogenerated
              "track" and "group" glyphs

  -subpart_callbacks Whether to invoke     false
              callbacks for subparts of
              the glyph.

  -box_subparts Return boxes around feature          false
               subparts rather than around the
               feature itself.

  -link, -title, -target
               These options are used when creating imagemaps
               for display on the web.  See L</"Creating Imagemaps">.

  -filter      Select which features to
               display. Must be a CODE reference.

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
value of 1 will "magically" look for tags of type "note" or
"description" and draw them if found, otherwise the source tag, if
any, will be displayed.  A code reference will be invoked to calculate
the description on the fly.  Anything else will be treated as a string
and used verbatim.

B<Connectors:> A glyph can contain subglyphs, recursively.  The top
level glyph is the track, which contains one or more groups, which
contain features, which contain subfeatures, and so forth.  By
default, the "group" glyph draws dotted lines between each of its
subglyphs, the "segment" glyph draws a solid line between each of its
subglyphs, and the "transcript" and "transcript2" glyphs draw
hat-shaped lines between their subglyphs.  All other glyphs do not
connect their components.  You can override this behavior by providing 
a -connector option, to explicitly set the type of connector.  Valid
options are:


   "hat"     an upward-angling conector
   "solid"   a straight horizontal connector
   "quill"   a decorated line with small arrows indicating strandedness
             (like the UCSC Genome Browser uses)
   "dashed"  a horizontal dashed line.

The B<-connector_color> option controls the color of the connector, if
any.

B<Collision control:> The B<-bump> argument controls what happens when
glyphs collide.  By default, they will simply overlap (value 0).  A
-bump value of +1 will cause overlapping glyphs to bump downwards
until there is room for them.  A -bump value of -1 will cause
overlapping glyphs to bump upwards.  You may also provide a -bump
value of +2 or -2 to activate a very simple type of collision control
in which each feature occupies its own line.  This is useful for
showing dense, nearly-full length features such as similarity hits.
The bump argument can also be a code reference; see below.

If you would like to see more horizontal whitespace between features
that occupy the same line, you can specify it with the B<-hbumppad>
option.  Positive values increase the amount of whitespace between
features.  Negative values decrease the whitespace.

B<Keys:> The -key argument declares that the track is to be shown in a
key appended to the bottom of the image.  The key contains a picture
of a glyph and a label describing what the glyph means.  The label is
specified in the argument to -key.

B<box_subparts:> Ordinarily, when you invoke the boxes() methods to
retrieve the rectangles surrounding the glyphs (which you need to do
to create clickable imagemaps, for example), the rectangles will
surround the top level features.  If you wish for the rectangles to
surround subpieces of the glyph, such as the exons in a transcript,
set box_subparts to a true numeric value. The value you specify will
control the number of levels of subfeatures that the boxes will
descend into. For example, if using the "gene" glyph, set
-box_subparts to 2 to create boxes for the whole gene (level 0), the
mRNAs (level 1) and the exons (level 2).

B<part_labels:> If set to true, each subpart of a multipart feature
will be labeled with a number starting with 1 at the 5'-most
part. This is useful for counting exons. You can pass a callback to
this argument; the part number and the total number of parts will be
arguments three and four. For example, to label the exons as "exon 1",
"exon 2" and so on:

 -part_labels  =>  sub {
		     my ($feature,undef,$partno) = @_;
		     return 'exon '.($partno+1);
	           }

The B<-label> argument must also be true.

B<part_labels_merge:> If true, changes the behavior of -part_labels so
that features that abut each other without a gap are treated as a
single feature. Useful if you want to count the UTR and CDS segments
of an exon as a single unit, and the default for transcript glyphs.

B<strand_arrow:> If set to true, some glyphs will indicate their
strandedness, usually by drawing an arrow.  For this to work, the
Bio::SeqFeature must have a strand of +1 or -1.  The glyph will ignore
this directive if the underlying feature has a strand of zero or
undef.

B<sort_order>: By default, features are drawn with a layout based only on the
position of the feature, assuring a maximal "packing" of the glyphs
when bumped.  In some cases, however, it makes sense to display the
glyphs sorted by score or some other comparison, e.g. such that more
"important" features are nearer the top of the display, stacked above
less important features.  The -sort_order option allows a few
different built-in values for changing the default sort order (which
is by "left" position): "low_score" (or "high_score") will cause
features to be sorted from lowest to highest score (or vice versa).
"left" (or "default") and "right" values will cause features to be
sorted by their position in the sequence.  "longest" (or "shortest")
will cause the longest (or shortest) features to be sorted first, and
"strand" will cause the features to be sorted by strand: "+1"
(forward) then "0" (unknown, or NA) then "-1" (reverse).

In all cases, the "left" position will be used to break any ties.  To
break ties using another field, options may be strung together using a
"|" character; e.g. "strand|low_score|right" would cause the features
to be sorted first by strand, then score (lowest to highest), then by
"right" position in the sequence.

Finally, a subroutine coderef with a $$ prototype can be provided.  It
will receive two B<glyph> as arguments and should return -1, 0 or 1
(see Perl's sort() function for more information).  For example, to
sort a set of database search hits by bits (stored in the features'
"score" fields), scaled by the log of the alignment length (with
"start" position breaking any ties):

  sort_order = sub ($$) {
    my ($glyph1,$glyph2) = @_;
    my $a = $glyph1->feature;
    my $b = $glyph2->feature;
    ( $b->score/log($b->length)
          <=>
      $a->score/log($a->length) )
          ||
    ( $a->start <=> $b->start )
  }

It is important to remember to use the $$ prototype as shown in the
example.  Otherwise Bio::Graphics will quit with an exception. The
arguments are subclasses of Bio::Graphics::Glyph, not the features
themselves.  While glyphs implement some, but not all, of the feature
methods, to be safe call the two glyphs' feature() methods in order to
convert them into the actual features.

The '-always_sort' option, if true, will sort features even if bumping
is turned off.  This is useful if you would like overlapping features
to stack in a particular order.  Features towards the end of the list
will overlay those towards the beginning of the sort order.

B<bump_limit>: When bumping is chosen, colliding features will
ordinarily move upward or downward without limit.  When many features
collide, this can lead to excessively high images.  You can limit the
number of levels that features will bump by providing a numeric
B<bump_limit> option.

The B<-filter> option, which must be a CODE reference, will be invoked
once for each feature prior to rendering it. The coderef will receive
the feature as its single option and should return true if the feature
is to be shown and false otherwise.

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
object itself.  The latter fields are useful in the case of treating
the first or last subfeature differently, such as using a different
color for the terminal exon of a gene.  Usually you will only need to
examine the first argument.  This example shows a callback examining
the score() attribute of a feature (possibly a BLAST hit) and return
the color "red" for high-scoring features, and "green" for low-scoring
features:

  sub callback {
     my $feature = shift;
     if ($feature->score > 90) {
       return 'red';
     else {
       return 'green';
    }
  }

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

In particular, this means that in order to control the -bump option
with a callback, you should specify -all_callbacks=E<gt>1, and turn on
bumping when the callback is in the track or group glyphs.

The -subpart_callbacks options is similar, except that when this is
set to true callbacks are invoked for the main glyph and its
subparts. This option only affects the -label and -description
options.

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
   flip()             Get/set coordinate flipping
   pad_top()	      Get/set top padding
   pad_left()	      Get/set left padding
   pad_bottom()	      Get/set bottom padding
   pad_right()	      Get/set right padding
   start()            Get the start of the sequence (bp; read only)
   end()              Get the end of the sequence (bp; read only)
   left()             Get the left side of the drawing area (pixels; read only)
   right()            Get the right side of the drawing area (pixels; read only)

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

=head2 Creating Imagemaps

You may wish to use Bio::Graphics to create clickable imagemaps for
display on the web.  The main method for achieving this is
image_and_map().  Under special circumstances you may instead wish to
call either or both of create_web_image() and create_web_map().

Here is a synopsis of how to use image_and_map() in a CGI script,
using CGI.pm calls to provide the HTML scaffolding:

   print h2('My Genome');

   my ($url,$map,$mapname) =
       $panel->image_and_map(-root => '/var/www/html',
                             -url  => '/tmpimages',
                             -link => 'http://www.google.com/search?q=$name');

   print img({-src=>$url,-usemap=>"#$mapname"});

   print $map;

We call image_and_map() with various arguments (described below) to
generate a three element list consisting of the URL at which the image
can be accessed, an HTML fragment containing the clickable imagemap
data, and the name of the map.  We print out an E<lt>imageE<gt> tag
that uses the URL of the map as its src attribute and the name of the
map as the value of its usemap attribute.  It is important to note
that we must put a "#" in front of the name of the map in order to
indicate that the map can be found in the same document as the
E<lt>imageE<gt> tag.  Lastly, we print out the map itself.

=over 4

=item ($url,$map,$mapname) = $panel-E<gt>image_and_map(@options)

Create the image in a web-accessible directory and return its URL, its
clickable imagemap, and the name of the imagemap.  The following
options are recognized:

 Option        Description
 ------        -----------

 -url          The URL to store the image at.


 -root         The directory path that should be appended to the
               start of -url in order to obtain a physical
               directory path.
 -link         A string pattern or coderef that will be used to
               generate the outgoing hypertext links for the imagemap.

 -title        A string pattern or coderef that will be used to
               generate the "title" tags of each element in the imagemap
               (these appear as popup hint boxes in certain browsers).

 -target       A string pattern or coderef that will be used to
               generate the window target for each element.  This can
               be used to pop up a new window when the user clicks on
               an element.

 -mapname      The name to use for the E<lt>mapE<gt> tag.  If not provided,
               a unique one will be autogenerated for you.

This method returns a three element list consisting of the URL at
which the image has been written to, the imagemap HTML, and the name
of the map.  Usually you will incorporate this information into an
HTML document like so:

  my ($url,$map,$mapname) =
          $panel->image_and_map(-link=>'http://www.google.com/searche?q=$name');
  print qq(<img src="$url" usemap="#$map">),"\n";
  print $map,"\n";

=item $url = $panel-E<gt>create_web_image($url,$root)

Create the image, write it into the directory indicated by
concatenating $root and $url (i.e. "$root/$url"), and return $url.

=item $map = $panel-E<gt>create_web_map('mapname',$linkrule,$titlerule,$targetrule)

Create a clickable imagemap named "mapname" using the indicated rules
to generate the hypertext links, the element titles, and the window
targets for the graphical elements.  Return the HTML for the map,
including the enclosing E<lt>mapE<gt> tag itself.

=back

To use this method effectively, you will need a web server and an
image directory in the document tree that is writable by the web
server user.  For example, if your web server's document root is
located at /var/www/html, you might want to create a directory named
"tmpimages" for this purpose:

  mkdir /var/www/html/tmpimages
  chmod 1777 /var/www/html/tmpimages

The 1777 privilege will allow anyone to create files and
subdirectories in this directory, but only the owner of the file will
be able to delete it.

When you call image_and_map(), you must provide it with two vital
pieces of information: the URL of the image directory and the physical
location of the web server's document tree.  In our example, you would
call:

  $panel->image_and_map(-root => '/var/www/html',-url=>'/tmpimages');

If you are working with virtual hosts, you might wish to provide the
hostname:portnumber part of the URL.  This will work just as well:

  $panel->image_and_map(-root => '/var/www/html',
                        -url  => 'http://myhost.com:8080/tmpimages');

If you do not provide the -root argument, the method will try to
figure it out from the DOCUMENT_ROOT environment variable.  If you do
not provide the -url argument, the method will assume "/tmp".

During execution, the image_and_map() method will generate a unique
name for the image using the Digest::MD5 module.  You can get this
module on CPAN and it B<must> be installed in order to use
image_and_map().  The imagename will be a long hexadecimal string such
as "e7457643f12d413f20843d4030c197c6.png".  Its URL will be
/tmpimages/e7457643f12d413f20843d4030c197c6.png, and its physical path
will be /var/www/html/tmpimages/e7457643f12d413f20843d4030c197c6.png

In addition to providing directory information, you must also tell
image_and_map() how to create outgoing links for each graphical
feature, and, optionally, how to create the "hover title" (the popup
yellow box displayed by most modern browsers), and the name of the
window or frame to link to when the user clicks on it.

There are three ways to specify the link destination:

=over 4

=item 1.

By configuring one or more tracks with a -link argument.

=item 2.

By configuring the panel with a -link argument.

=item 3.

By passing a -link argument in the call to image_and_map().

=back

The -link argument can be either a string or a coderef.  If you pass a
string, it will be interpreted as a URL pattern containing runtime
variables.  These variables begin with a dollar sign ($), and are
replaced at run time with the information relating to the selected
annotation.  Recognized variables include:

     $name        The feature's name (display name)
     $id          The feature's id (eg, PK from a database)
     $class       The feature's class (group class)
     $method      The feature's method (same as primary tag)
     $source      The feature's source
     $ref         The name of the sequence segment (chromosome, contig)
                     on which this feature is located
     $description The feature's description (notes)
     $start       The start position of this feature, relative to $ref
     $end         The end position of this feature, relative to $ref
     $segstart    The left end of $ref displayed in the detailed view
     $segend      The right end of $ref displayed in the detailed view

For example, to link each feature to a Google search on the feature's
description, use the argument:

  -link => 'http://www.google.com/search?q=$description'

Be sure to use single quotes around the pattern, or Perl will attempt
to perform variable interpretation before image_and_map() has a chance
to work on it.

You may also pass a code reference to -link, in which case the code
will be called every time a URL needs to be generated for the
imagemap.  The subroutine will be called with two arguments, the
feature and the Bio::Graphics::Panel object, and it should return the
URL to link to, or an empty string if a link is not desired. Here is a
simple example:

  -link => sub {
         my ($feature,$panel) = @_;
         my $type = $feature->primary_tag;
         my $name = $feature->display_name;
         if ($primary_tag eq 'clone') {
            return "http://www.google.com/search?q=$name";
         } else {
            return "http://www.yahoo.com/search?p=$name";
         }

The -link argument cascades. image_and_map() will first look for a
-link option in the track configuration, and if that's not found, it
will look in the Panel configuration (created during
Bio::Graphics::Panel-E<gt>new). If no -link configuration option is found
in either location, then image_and_map() will use the value of -link
passed in its argument list, if any.

The -title and -target options behave in a similar manner to -link.
-title is used to assign each feature "title" and "alt" attributes.
The "title" attribute is used by many browsers to create a popup hints
box when the mouse hovers over the feature's glyph for a preset length
of time, while the "alt" attribute is used to create navigable menu
items for the visually impaired.  As with -link, you can set the title
by passing either a substitution pattern or a code ref, and the -title
option can be set in the track, the panel, or the method call itself
in that order of priority.

If not provided, image_and_map() will autogenerate its own title in
the form "E<lt>methodE<gt> E<lt>display_nameE<gt> E<lt>seqidE<gt>:start..end".

The -target option can be used to specify the window or frame that
clicked features will link to.  By default, when the user clicks on a
feature, the loaded URL will replace the current page.  You can modify
this by providing -target with the name of a preexisting or new window
name in order to create effects like popup windows, multiple frames,
popunders and the like.  The value of -target follows the same rules
as -title and -link, including variable substitution and the use of
code refs.

NOTE: Each time you call image_and_map() it will generate a new image
file.  Images that are identical to an earlier one will reuse the same
name, but those that are different, even by one pixel, will result in
the generation of a new image.  If you have limited disk space, you
might wish to check the images directory periodically and remove those
that have not been accessed recently.  The following cron script will
remove image files that haven't been accessed in more than 20 days.

30 2 * * * find /var/www/html/tmpimages -type f -atime +20 -exec rm {} \;

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
L<Bio::Graphics::Glyph::redgreen_box>,
L<Bio::Graphics::Glyph::ruler_arrow>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::translation>,
L<Bio::Graphics::Glyph::triangle>,
L<Bio::Graphics::Glyph::xyplot>,
L<Bio::Graphics::Glyph::whiskerplot>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>
L<GD::SVG>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

