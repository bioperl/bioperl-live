package Bio::Graphics::Glyph;
use GD;

use strict;
use Carp 'croak';
use constant BUMP_SPACING => 2; # vertical distance between bumped glyphs


my %LAYOUT_COUNT;

# the CM1 and CM2 constants control the size of the hash used to
# detect collisions.
use constant CM1 => 200; # big bin, x axis
use constant CM2 => 50;  # big bin, y axis
use constant CM3 => 50;  # small bin, x axis
use constant CM4 => 50;  # small bin, y axis

use constant QUILL_INTERVAL => 8;  # number of pixels between Jim Kent style intron "quills"

# a bumpable graphical object that has bumpable graphical subparts

# args:  -feature => $feature_object (may contain subsequences)
#        -factory => $factory_object (called to create glyphs for subsequences)
# In this scheme, the factory decides based on stylesheet information what glyph to
# draw and what configurations options to us. This allows for heterogeneous tracks.
sub new {
  my $class = shift;
  my %arg = @_;

  my $feature = $arg{-feature} or die "No feature";
  my $factory = $arg{-factory} || $class->default_factory;
  my $level   = $arg{-level} || 0;
  my $flip    = $arg{-flip};

  my $self = bless {},$class;
  $self->{feature} = $feature;
  $self->{factory} = $factory;
  $self->{level}   = $level;
  $self->{flip}++  if $flip;
  $self->{top} = 0;

  my @subglyphs;
  my @subfeatures = $self->subseq($feature);

  if (@subfeatures) {

    # dynamic glyph resolution
    @subglyphs = map { $_->[0] }
          sort { $a->[1] <=> $b->[1] }
             map { [$_, $_->left ] } 
    $factory->make_glyph($level+1,@subfeatures);

    $self->{parts}   = \@subglyphs;
  }

  my ($start,$stop) = ($self->start, $self->stop);
  if (defined $start && defined $stop) {
    ($start,$stop) = ($stop,$start) if $start > $stop;  # sheer paranoia
    # the +1 here is critical for allowing features to meet nicely at nucleotide resolution
    my ($left,$right) = $factory->map_pt($start,$stop+1);
    $self->{left}    = $left;
    $self->{width}   = $right - $left + 1;
  }
  if (@subglyphs) {
      my $l            = $subglyphs[0]->left;
      $self->{left}    = $l if !defined($self->{left}) || $l < $self->{left};
      my $right        = (
			  sort { $b<=>$a } 
			  map {$_->right} @subglyphs)[0];
      my $w            = $right - $self->{left} + 1;
      $self->{width}   = $w if !defined($self->{width}) || $w > $self->{width};
  }

  $self->{point} = $arg{-point} ? $self->height : undef;
  #Handle glyphs that don't actually fill their space, but merely mark a point.
  #They need to have their collision bounds altered.  We will (for now)
  #hard code them to be in the center of their feature.
# note: this didn't actually seem to work properly, all features were aligned on
# their right edges.  It works to do it in individual point-like glyphs such as triangle.
#  if($self->option('point')){
#    my ($left,$right) = $factory->map_pt($self->start,$self->stop);
#    my $center = int(($left+$right)/2 + 0.5);

#    $self->{width} = $self->height;
#    $self->{left}  = $center - ($self->{width});
#    $self->{right} = $center + ($self->{width});
#  }

  return $self;
}

sub parts      {
  my $self = shift;
  return unless $self->{parts};
  return wantarray ? @{$self->{parts}} : $self->{parts};
}

sub feature { shift->{feature} }
sub factory { shift->{factory} }
sub panel   { shift->factory->panel }
sub point   { shift->{point}   }
sub scale   { shift->factory->scale }
sub start   {
  my $self = shift;
  return $self->{start} if exists $self->{start};
  $self->{start} = $self->{flip} ? $self->panel->end + 1 - $self->{feature}->end : $self->{feature}->start;

  # handle the case of features whose endpoints are undef
  # (this happens with wormbase clones where one or more clone end is not defined)
  # in this case, we set the start to one minus the beginning of the panel
  $self->{start} = $self->panel->offset - 1 unless defined $self->{start};

  return $self->{start};
}
sub stop    {
  my $self = shift;
  return $self->{stop} if exists $self->{stop};
  $self->{stop} = $self->{flip} ?  $self->panel->end + 1 - $self->{feature}->start : $self->{feature}->end;

  # handle the case of features whose endpoints are undef
  # (this happens with wormbase clones where one or more clone end is not defined)
  # in this case, we set the start to one plus the end of the panel
  $self->{stop} = $self->panel->offset + $self->panel->length + 1 unless defined $self->{stop};

  return $self->{stop}
}
sub end     { shift->stop }
sub length { my $self = shift; $self->stop - $self->start };
sub score {
    my $self = shift;
    return $self->{score} if exists $self->{score};
    return $self->{score} = ($self->{feature}->score || 0);
}
sub strand {
    my $self = shift;
    return $self->{strand} if exists $self->{strand};
    return $self->{strand} = ($self->{feature}->strand || 0);
}
sub map_pt  { shift->{factory}->map_pt(@_) }
sub map_no_trunc { shift->{factory}->map_no_trunc(@_) }

# add a feature (or array ref of features) to the list
sub add_feature {
  my $self       = shift;
  my $factory    = $self->factory;
  for my $feature (@_) {
    if (ref $feature eq 'ARRAY') {
      $self->add_group(@$feature);
    } else {
      push @{$self->{parts}},$factory->make_glyph(0,$feature);
    }
  }
}

# link a set of features together so that they bump as a group
sub add_group {
  my $self = shift;
  my @features = ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;
  my $f    = Bio::Graphics::Feature->new(
					 -segments=>\@features,
					 -type => 'group'
					);
  $self->add_feature($f);
}

sub top {
  my $self = shift;
  my $g = $self->{top};
  $self->{top} = shift if @_;
  $g;
}
sub left {
  my $self = shift;
  return $self->{left} - $self->pad_left;
}
sub right {
  my $self = shift;
  return $self->left + $self->layout_width - 1;
}
sub bottom {
  my $self = shift;
  $self->top + $self->layout_height - 1;
}
sub height {
  my $self = shift;
  return $self->{height} if exists $self->{height};
  my $baseheight = $self->option('height');  # what the factory says
  return $self->{height} = $baseheight;
}
sub width {
  my $self = shift;
  my $g = $self->{width};
  $self->{width} = shift if @_;
  $g;
}
sub layout_height {
  my $self = shift;
  return $self->layout;
}
sub layout_width {
  my $self = shift;
  return $self->width + $self->pad_left + $self->pad_right;
}

# returns the rectangle that surrounds the physical part of the
# glyph, excluding labels and other "extra" stuff
sub calculate_boundaries {return shift->bounds(@_);}

sub bounds {
  my $self = shift;
  my ($dx,$dy) = @_;
  $dx += 0; $dy += 0;
  ($dx + $self->{left},
   $dy + $self->top    + $self->pad_top,
   $dx + $self->{left} + $self->{width} - 1,
   $dy + $self->bottom - $self->pad_bottom);
}


sub box {
  my $self = shift;
  return ($self->left,$self->top,$self->right,$self->bottom);
}


sub unfilled_box {
  my $self = shift;
  my $gd   = shift;
  my ($x1,$y1,$x2,$y2,$fg,$bg) = @_;

  my $linewidth = $self->option('linewidth') || 1;

  unless ($fg) {
      $fg ||= $self->fgcolor;
  $fg = $self->set_pen($linewidth,$fg) if $linewidth > 1;
  }

  unless ($bg) {
      $bg ||= $self->bgcolor;
      $bg = $self->set_pen($linewidth,$bg) if $linewidth > 1;
  }

  # draw a box
  $gd->rectangle($x1,$y1,$x2,$y2,$fg);

  # if the left end is off the end, then cover over
  # the leftmost line
  my ($width) = $gd->getBounds;

  $gd->line($x1,$y1+$linewidth,$x1,$y2-$linewidth,$bg)
    if $x1 < $self->panel->pad_left;

  $gd->line($x2,$y1+$linewidth,$x2,$y2-$linewidth,$bg)
    if $x2 > $width - $self->panel->pad_right;
}


# return boxes surrounding each part
sub boxes {
  my $self = shift;
  my ($left,$top) = @_;
  $top  += 0; $left += 0;
  my @result;

  $self->layout;
  my @parts = $self->parts;
  @parts    = $self if !@parts && $self->option('box_subparts') && $self->level>0;

  for my $part ($self->parts) {
    if (eval{$part->feature->primary_tag} eq 'group' or
	($part->level == 0 && $self->option('box_subparts'))) {
      push @result,$part->boxes($left+$self->left+$self->pad_left,$top+$self->top+$self->pad_top);
    } else {
      my ($x1,$y1,$x2,$y2) = $part->box;
      push @result,[$part->feature,$x1,$top+$self->top+$self->pad_top+$y1,
		                   $x2,$top+$self->top+$self->pad_top+$y2];
    }
  }
  return wantarray ? @result : \@result;
}

# this should be overridden for labels, etc.
# allows glyph to make itself thicker or thinner depending on
# domain-specific knowledge
sub pad_top {
  my $self = shift;
  return 0;
}
sub pad_bottom {
  my $self = shift;
  return 0;
}
sub pad_left {
  my $self = shift;
  return 0;
}
sub pad_right {
  my $self = shift;
# this shouldn't be necessary
  my @parts = $self->parts or return 0;
  my $max = 0;
  foreach (@parts) {
    my $pr = $_->pad_right;
    $max = $pr if $max < $pr;
  }
  $max;
}

# move relative to parent
sub move {
  my $self = shift;
  my ($dx,$dy) = @_;
  $self->{left} += $dx;
  $self->{top}  += $dy;

  # because the feature parts use *absolute* not relative addressing
  # we need to move each of the parts horizontally, but not vertically
  $_->move($dx,0) foreach $self->parts;
}

# get an option
sub option {
  my $self = shift;
  my $option_name = shift;
  my $factory = $self->factory;
  return unless $factory;
  $factory->option($self,$option_name,@{$self}{qw(partno total_parts)});
}

# set an option globally
sub configure {
  my $self = shift;
  my $factory = $self->factory;
  my $option_map = $factory->option_map;
  while (@_) {
    my $option_name  = shift;
    my $option_value = shift;
    ($option_name = lc $option_name) =~ s/^-//;
    $option_map->{$option_name} = $option_value;
  }
}

# some common options
sub color {
  my $self = shift;
  my $color = shift;
  my $index = $self->option($color);
  # turn into a color index
  return $self->factory->translate_color($index) if defined $index;
  return 0;
}

sub connector {
  return shift->option('connector',@_);
}

# return value:
#              0    no bumping
#              +1   bump down
#              -1   bump up
sub bump {
  my $self = shift;
  return $self->option('bump');
}

# we also look for the "color" option for Ace::Graphics compatibility
sub fgcolor {
  my $self = shift;
  my $color = $self->option('fgcolor');
  my $index = defined $color ? $color : $self->option('color');
  $index = 'black' unless defined $index;
  $self->factory->translate_color($index);
}

#add for compatibility
sub fillcolor {
    my $self = shift;
    return $self->bgcolor;
}

# we also look for the "background-color" option for Ace::Graphics compatibility
sub bgcolor {
  my $self = shift;
  my $bgcolor = $self->option('bgcolor');
  my $index = defined $bgcolor ? $bgcolor : $self->option('fillcolor');
  $index = 'white' unless defined $index;
  $self->factory->translate_color($index);
}
sub font {
  my $self = shift;
  my $font = $self->option('font');
  unless (UNIVERSAL::isa($font,'GD::Font')) {
    my $ref    = {
		  gdTinyFont  => gdTinyFont,
		  gdSmallFont => gdSmallFont,
		  gdMediumBoldFont => gdMediumBoldFont,
		  gdLargeFont => gdLargeFont,
		  gdGiantFont => gdGiantFont};
    my $gdfont = $ref->{$font} || $font;
    $self->configure(font=>$gdfont);
    return $gdfont;
  }
  return $font;
}
sub fontcolor {
  my $self = shift;
  my $fontcolor = $self->color('fontcolor');
  return defined $fontcolor ? $fontcolor : $self->fgcolor;
}
sub font2color {
  my $self = shift;
  my $font2color = $self->color('font2color');
  return defined $font2color ? $font2color : $self->fgcolor;
}
sub tkcolor { # "track color"
  my $self = shift;
  $self->option('tkcolor') or return;
  return $self->color('tkcolor')
}
sub connector_color {
  my $self = shift;
  $self->color('connector_color') || $self->fgcolor;
}

sub layout_sort {

    my $self = shift;
    my $sortfunc;

    my $opt = $self->option("sort_order");
    if (!$opt) {
       $sortfunc = eval 'sub { $a->left <=> $b->left }';
    } elsif (ref $opt eq 'CODE') {
       $sortfunc = $opt;
    } elsif ($opt =~ /^sub\s+\{/o) {
       $sortfunc = eval $opt;
    } else {
       # build $sortfunc for ourselves:
       my @sortbys = split(/\s*\|\s*/o, $opt);
       $sortfunc = 'sub { ';
       my $sawleft = 0;

       # not sure I can make this schwartzian transfored
       for my $sortby (@sortbys) {
	 if ($sortby eq "left" || $sortby eq "default") {
	   $sortfunc .= '($a->left <=> $b->left) || ';
	   $sawleft++;
	 } elsif ($sortby eq "right") {
	   $sortfunc .= '($a->right <=> $b->right) || ';
	 } elsif ($sortby eq "low_score") {
	   $sortfunc .= '($a->score <=> $b->score) || ';
	 } elsif ($sortby eq "high_score") {
	   $sortfunc .= '($b->score <=> $a->score) || ';
	 } elsif ($sortby eq "longest") {
	   $sortfunc .= '(($b->length) <=> ($a->length)) || ';
	 } elsif ($sortby eq "shortest") {
	   $sortfunc .= '(($a->length) <=> ($b->length)) || ';
	 } elsif ($sortby eq "strand") {
	   $sortfunc .= '($b->strand <=> $a->strand) || ';
	 } elsif ($sortby eq "name") {
	   $sortfunc .= '($a->feature->display_name cmp $b->feature->display_name) || ';
	 }
       }
       unless ($sawleft) {
           $sortfunc .= ' ($a->left <=> $b->left) ';
       } else {
           $sortfunc .= ' 0';
       }
       $sortfunc .= '}';
       $sortfunc = eval $sortfunc;
    }

    # cache this
    # $self->factory->set_option(sort_order => $sortfunc);

    return sort $sortfunc @_;
}

# handle collision detection
sub layout {
  my $self = shift;
  return $self->{layout_height} if exists $self->{layout_height};

  my @parts = $self->parts;
  return $self->{layout_height}
    = $self->height + $self->pad_top + $self->pad_bottom unless @parts;

  my $bump_direction = $self->bump;
  my $bump_limit = $self->option('bump_limit') || -1;

  $_->layout foreach @parts;  # recursively lay out

  # no bumping requested, or only one part here
  if (@parts == 1 || !$bump_direction) {
    my $highest = 0;
    foreach (@parts) {
      my $height = $_->layout_height;
      $highest   = $height > $highest ? $height : $highest;
    }
    return $self->{layout_height} = $highest + $self->pad_top + $self->pad_bottom;
  }

  my (%bin1,%bin2);
  for my $g ($self->layout_sort(@parts)) {

    my $pos = 0;
    my $bumplevel = 0;
    my $left   = $g->left;
    my $right  = $g->right;
    my $height = $g->{layout_height};

    while (1) {

      # stop bumping if we've gone too far down
      if ($bump_limit > 0 && $bumplevel++ >= $bump_limit) {
	$g->{overbumped}++;  # this flag can be used to suppress label and description
	foreach ($g->parts) {
	  $_->{overbumped}++;
	}
	last;
      }

      # look for collisions
      my $bottom = $pos + $height;
      $self->collides(\%bin1,CM1,CM2,$left,$pos,$right,$bottom) or last;
      my $collision = $self->collides(\%bin2,CM3,CM4,$left,$pos,$right,$bottom) or last;

      if ($bump_direction > 0) {
	$pos += $collision->[3]-$collision->[1] + BUMP_SPACING;    # collision, so bump

      } else {
	$pos -= BUMP_SPACING;
      }

    }

    $g->move(0,$pos);
    $self->add_collision(\%bin1,CM1,CM2,$left,$g->top,$right,$g->bottom);
    $self->add_collision(\%bin2,CM3,CM4,$left,$g->top,$right,$g->bottom);
  }

  # If -1 bumping was allowed, then normalize so that the top glyph is at zero
  if ($bump_direction < 0) {
    my $topmost;
    foreach (@parts) {
      my $top  = $_->top;
      $topmost = $top if !defined($topmost) or $top < $topmost;
    }
    my $offset = - $topmost;
    $_->move(0,$offset) foreach @parts;
  }

  # find new height
  my $bottom = 0;
  foreach (@parts) {
    $bottom = $_->bottom if $_->bottom > $bottom;
  }
  return $self->{layout_height} = $self->pad_bottom + $self->pad_top + $bottom - $self->top  + 1;
}

# the $%occupied structure is a hash of {left,top} = [left,top,right,bottom]
sub collides {
  my $self = shift;
  my ($occupied,$cm1,$cm2,$left,$top,$right,$bottom) = @_;
  my @keys = $self->_collision_keys($cm1,$cm2,$left,$top,$right,$bottom);
  my $collides = 0;
  for my $k (@keys) {
    next unless exists $occupied->{$k};
    for my $bounds (@{$occupied->{$k}}) {
      my ($l,$t,$r,$b) = @$bounds;
      next unless $right >= $l and $left <= $r and $bottom >= $t and $top <= $b;
      $collides = $bounds;
      last;
    }
  }
  $collides;
}

sub add_collision {
  my $self = shift;
  my ($occupied,$cm1,$cm2,$left,$top,$right,$bottom) = @_;
  my $value = [$left,$top,$right+2,$bottom];
  my @keys = $self->_collision_keys($cm1,$cm2,@$value);
  push @{$occupied->{$_}},$value foreach @keys;
}

sub _collision_keys {
  my $self = shift;
  my ($binx,$biny,$left,$top,$right,$bottom) = @_;
  my @keys;
  my $bin_left   = int($left/$binx);
  my $bin_right  = int($right/$binx);
  my $bin_top    = int($top/$biny);
  my $bin_bottom = int($bottom/$biny);
  for (my $x=$bin_left;$x<=$bin_right; $x++) {
    for (my $y=$bin_top;$y<=$bin_bottom; $y++) {
      push @keys,join(',',$x,$y);
    }
  }
  @keys;
}

sub draw {
  my $self = shift;
  my $gd = shift;
  my ($left,$top,$partno,$total_parts) = @_;

  local($self->{partno},$self->{total_parts});
  @{$self}{qw(partno total_parts)} = ($partno,$total_parts);

  my $connector =  $self->connector;
  if (my @parts = $self->parts) {

    # invoke sorter if use wants to sort always and we haven't already sorted
    # during bumping.
    @parts = $self->layout_sort(@parts) if !$self->bump && $self->option('always_sort');

    my $x = $left;
    my $y = $top  + $self->top + $self->pad_top;
    $self->draw_connectors($gd,$x,$y) if $connector && $connector ne 'none';

    my $last_x;
    for (my $i=0; $i<@parts; $i++) {
      # lie just a little bit to avoid lines overlapping and
      # make the picture prettier
      my $fake_x = $x;
      $fake_x-- if defined $last_x && $parts[$i]->left - $last_x == 1;
      $parts[$i]->draw($gd,$fake_x,$y,$i,scalar(@parts));
      $last_x = $parts[$i]->right;
    }
  }

  else {  # no part
    $self->draw_connectors($gd,$left,$top)
      if $connector && $connector ne 'none' && $self->{level} == 0;
    $self->draw_component($gd,$left,$top);
  }
}

# the "level" is the level of testing of the glyph
# groups are level -1, top level glyphs are level 0, subcomponents are level 1 and so forth.
sub level {
  shift->{level};
}

sub draw_connectors {
  my $self = shift;
  return if $self->{overbumped};
  my $gd = shift;
  my ($dx,$dy) = @_;
  my @parts = sort { $a->left <=> $b->left } $self->parts;
  for (my $i = 0; $i < @parts-1; $i++) {
    $self->_connector($gd,$dx,$dy,$parts[$i]->bounds,$parts[$i+1]->bounds);
  }

  # extra connectors going off ends
  if (@parts) {
    my($x1,$y1,$x2,$y2) = $self->bounds(0,0);
    my($xl,$xt,$xr,$xb) = $parts[0]->bounds;
    $self->_connector($gd,$dx,$dy,$x1,$xt,$x1,$xb,$xl,$xt,$xr,$xb)      if $x1 < $xl;
    my ($xl2,$xt2,$xr2,$xb2) = $parts[-1]->bounds;
    $self->_connector($gd,$dx,$dy,$parts[-1]->bounds,$x2,$xt2,$x2,$xb2) if $x2 > $xr;
  }

}

sub _connector {
  my $self = shift;
  my ($gd,
      $dx,$dy,
      $xl,$xt,$xr,$xb,
      $yl,$yt,$yr,$yb) = @_;
  my $left   = $dx + $xr;
  my $right  = $dx + $yl;
  my $top1     = $dy + $xt;
  my $bottom1  = $dy + $xb;
  my $top2     = $dy + $yt;
  my $bottom2  = $dy + $yb;
  # restore this comment if you don't like the group dash working
  # its way backwards.
  return if $right-$left < 1 && !$self->isa('Bio::Graphics::Glyph::group');

  $self->draw_connector($gd,
			$top1,$bottom1,$left,
			$top2,$bottom2,$right,
		       );
}

sub draw_connector {
  my $self   = shift;
  my $gd     = shift;

  my $color          = $self->connector_color;
  my $connector_type = $self->connector or return;

  if ($connector_type eq 'hat') {
    $self->draw_hat_connector($gd,$color,@_);
  } elsif ($connector_type eq 'solid') {
    $self->draw_solid_connector($gd,$color,@_);
  } elsif ($connector_type eq 'dashed') {
    $self->draw_dashed_connector($gd,$color,@_);
  } elsif ($connector_type eq 'quill') {
    $self->draw_quill_connector($gd,$color,@_);
  } else {
    ; # draw nothing
  }
}

sub draw_hat_connector {
  my $self = shift;
  my $gd   = shift;
  my $color = shift;
  my ($top1,$bottom1,$left,$top2,$bottom2,$right) = @_;

  my $center1  = ($top1 + $bottom1)/2;
  my $quarter1 = $top1 + ($bottom1-$top1)/4;
  my $center2  = ($top2 + $bottom2)/2;
  my $quarter2 = $top2 + ($bottom2-$top2)/4;

  if ($center1 != $center2) {
    $self->draw_solid_connector($gd,$color,@_);
    return;
  }

  if ($right - $left > 4) {  # room for the inverted "V"
      my $middle = $left + int(($right - $left)/2);
      $gd->line($left,$center1,$middle,$top1,$color);
      $gd->line($middle,$top1,$right-1,$center1,$color);
    } elsif ($right-$left > 1) { # no room, just connect
      $gd->line($left,$quarter1,$right-1,$quarter1,$color);
    }

}

sub draw_solid_connector {
  my $self = shift;
  my $gd   = shift;
  my $color = shift;
  my ($top1,$bottom1,$left,$top2,$bottom2,$right) = @_;

  my $center1  = ($top1 + $bottom1)/2;
  my $center2  = ($top2 + $bottom2)/2;

  $gd->line($left,$center1,$right,$center2,$color);
}

sub draw_dashed_connector {
  my $self = shift;
  my $gd   = shift;
  my $color = shift;
  my ($top1,$bottom1,$left,$top2,$bottom2,$right) = @_;

  my $center1  = ($top1 + $bottom1)/2;
  my $center2  = ($top2 + $bottom2)/2;

  $gd->setStyle($color,$color,gdTransparent,gdTransparent,);
  $gd->line($left,$center1,$right,$center2,gdStyled);
}

sub draw_quill_connector {
  my $self = shift;
  my $gd   = shift;
  my $color = shift;
  my ($top1,$bottom1,$left,$top2,$bottom2,$right) = @_;

  my $center1  = ($top1 + $bottom1)/2;
  my $center2  = ($top2 + $bottom2)/2;

  $gd->line($left,$center1,$right,$center2,$color);
  my $direction = $self->feature->strand;
  return unless $direction;

  if ($direction > 0) {
    my $start = $left+4;
    my $end   = $right-1;
    for (my $position=$start; $position <= $end; $position += QUILL_INTERVAL) {
      $gd->line($position,$center1,$position-2,$center1-2,$color);
      $gd->line($position,$center1,$position-2,$center1+2,$color);
    }
  } else {
    my $start = $left+1;
    my $end   = $right-4;
    for (my $position=$start; $position <= $end; $position += QUILL_INTERVAL) {
      $gd->line($position,$center1,$position+2,$center1-2,$color);
      $gd->line($position,$center1,$position+2,$center1+2,$color);
    }
  }
}

sub filled_box {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2,$bg,$fg) = @_;

  $bg ||= $self->bgcolor;
  $fg ||= $self->fgcolor;
  my $linewidth = $self->option('linewidth') || 1;

  $gd->filledRectangle($x1,$y1,$x2,$y2,$bg);

  $fg = $self->set_pen($linewidth,$fg) if $linewidth > 1;

  # draw a box
  $gd->rectangle($x1,$y1,$x2,$y2,$fg);

  # if the left end is off the end, then cover over
  # the leftmost line
  my ($width) = $gd->getBounds;

  $bg = $self->set_pen($linewidth,$bg) if $linewidth > 1;

  $gd->line($x1,$y1+$linewidth,$x1,$y2-$linewidth,$bg)
    if $x1 < $self->panel->pad_left;

  $gd->line($x2,$y1+$linewidth,$x2,$y2-$linewidth,$bg)
    if $x2 > $width - $self->panel->pad_right;
}

sub filled_oval {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2,$bg,$fg) = @_;
  my $cx = ($x1+$x2)/2;
  my $cy = ($y1+$y2)/2;

  $fg ||= $self->fgcolor;
  $bg ||= $self->bgcolor;
  my $linewidth = $self->linewidth;

  $fg = $self->set_pen($linewidth) if $linewidth > 1;
  $gd->arc($cx,$cy,$x2-$x1,$y2-$y1,0,360,$fg);

  # and fill it
  $gd->fill($cx,$cy,$bg);
}

sub oval {
  my $self = shift;
  my $gd = shift;
  my ($x1,$y1,$x2,$y2) = @_;
  my $cx = ($x1+$x2)/2;
  my $cy = ($y1+$y2)/2;

  my $fg = $self->fgcolor;
  my $linewidth = $self->linewidth;

  $fg = $self->set_pen($linewidth) if $linewidth > 1;
  $gd->arc($cx,$cy,$x2-$x1,$y2-$y1,0,360,$fg);
}

sub filled_arrow {
  my $self = shift;
  my $gd  = shift;
  my $orientation = shift;
  $orientation *= -1 if $self->{flip};

  my ($x1,$y1,$x2,$y2) = @_;

  my ($width) = $gd->getBounds;
  my $indent = $y2-$y1 < $x2-$x1 ? $y2-$y1 : ($x2-$x1)/2;

  return $self->filled_box($gd,@_)
    if ($orientation == 0)
      or ($x1 < 0 && $orientation < 0)
        or ($x2 > $width && $orientation > 0)
	  or ($indent <= 0)
	    or ($x2 - $x1 < 3);

  my $fg = $self->fgcolor;
  if ($orientation >= 0) {
    $gd->line($x1,$y1,$x2-$indent,$y1,$fg);
    $gd->line($x2-$indent,$y1,$x2,($y2+$y1)/2,$fg);
    $gd->line($x2,($y2+$y1)/2,$x2-$indent,$y2,$fg);
    $gd->line($x2-$indent,$y2,$x1,$y2,$fg);
    $gd->line($x1,$y2,$x1,$y1,$fg);
    my $left = $self->panel->left > $x1 ? $self->panel->left : $x1;
    $gd->fillToBorder($left+1,($y1+$y2)/2,$fg,$self->bgcolor);
  } else {
    $gd->line($x1,($y2+$y1)/2,$x1+$indent,$y1,$fg);
    $gd->line($x1+$indent,$y1,$x2,$y1,$fg);
    $gd->line($x2,$y2,$x1+$indent,$y2,$fg);
    $gd->line($x1+$indent,$y2,$x1,($y1+$y2)/2,$fg);
    $gd->line($x2,$y1,$x2,$y2,$fg);
    my $right = $self->panel->right < $x2 ? $self->panel->right : $x2;
    $gd->fillToBorder($right-1,($y1+$y2)/2,$fg,$self->bgcolor);
  }
}

sub linewidth {
  shift->option('linewidth') || 1;
}

sub fill {
  my $self = shift;
  my $gd   = shift;
  my ($x1,$y1,$x2,$y2) = @_;
  if ( ($x2-$x1) >= 2 && ($y2-$y1) >= 2 ) {
    $gd->fill($x1+1,$y1+1,$self->bgcolor);
  }
}
sub set_pen {
  my $self = shift;
  my ($linewidth,$color) = @_;
  $linewidth ||= $self->linewidth;
  $color     ||= $self->fgcolor;
  return $color unless $linewidth > 1;
  $self->panel->set_pen($linewidth,$color);
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my($x1,$y1,$x2,$y2) = $self->bounds(@_);

  # clipping
  my $panel = $self->panel;
  return unless $x2 >= $panel->left and $x1 <= $panel->right;

  if ($self->option('strand_arrow') || $self->option('stranded')) {
    $self->filled_arrow($gd,$self->feature->strand,
			$x1, $y1,
			$x2, $y2)
  } else {
    $self->filled_box($gd,
		      $x1, $y1,
		      $x2, $y2)
  }
}

# memoize _subseq -- it's a bottleneck with segments
sub subseq {
  my $self    = shift;
  my $feature = shift;
  return $self->_subseq($feature) unless ref $self;
  return @{$self->{cached_subseq}{$feature}} if $self->{cached_subseq}{$feature};
  my @ss = $self->_subseq($feature);
  $self->{cached_subseq}{$feature} = \@ss;
  @ss;
}

sub _subseq {
  my $class   = shift;
  my $feature = shift;
  return $feature->merged_segments         if $feature->can('merged_segments');
  return $feature->segments                if $feature->can('segments');
  my @split = eval { my $id   = $feature->location->seq_id;
		     my @subs = $feature->location->sub_Location;
		     grep {$id eq $_->seq_id} @subs};
  return @split if @split;
  return $feature->sub_SeqFeature          if $feature->can('sub_SeqFeature');
  return;
}

# synthesize a key glyph
sub keyglyph {
  my $self = shift;
  my $feature = $self->make_key_feature;
  my $factory = $self->factory->clone;
  $factory->set_option(label       => 1);
  $factory->set_option(description => 0);
  $factory->set_option(bump  => 0);
  $factory->set_option(connector  => 'solid');
  return $factory->make_glyph(0,$feature);
}

# synthesize a key glyph
sub make_key_feature {
  my $self = shift;

  my $scale = 1/$self->scale;  # base pairs/pixel

  # one segments, at pixels 0->80
  my $offset = $self->panel->offset;


  my $feature =
    Bio::Graphics::Feature->new(-start =>0 * $scale +$offset,
				-end   =>80*$scale+$offset,
				-name => $self->option('key'),
				-strand => '+1');
  return $feature;
}

sub all_callbacks {
  my $self = shift;
  my $track_level = $self->option('all_callbacks');
  return $track_level if defined $track_level;
  return $self->panel->all_callbacks;
}

sub default_factory {
  croak "no default factory implemented";
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph - Base class for Bio::Graphics::Glyph objects

=head1 SYNOPSIS

See L<Bio::Graphics::Panel>.

=head1 DESCRIPTION

Bio::Graphics::Glyph is the base class for all glyph objects.  Each
glyph is a wrapper around an Bio:SeqFeatureI object, knows how to
render itself on an Bio::Graphics::Panel, and has a variety of
configuration variables.

End developers will not ordinarily work directly with
Bio::Graphics::Glyph objects, but with Bio::Graphics::Glyph::generic
and its subclasses.  Similarly, most glyph developers will want to
subclass from Bio::Graphics::Glyph::generic because the latter
provides labeling and arrow-drawing facilities.

=head1 METHODS

This section describes the class and object methods for
Bio::Graphics::Glyph.

=head2 CONSTRUCTORS

Bio::Graphics::Glyph objects are constructed automatically by an
Bio::Graphics::Glyph::Factory, and are not usually created by
end-developer code.

=over 4

=item $glyph = Bio::Graphics::Glyph-E<gt>new(-feature=E<gt>$feature,-factory=E<gt>$factory)

Given a sequence feature, creates an Bio::Graphics::Glyph object to
display it.  The B<-feature> argument points to the Bio:SeqFeatureI
object to display, and B<-factory> indicates an
Bio::Graphics::Glyph::Factory object from which the glyph will fetch
all its run-time configuration information.  Factories are created and
manipulated by the Bio::Graphics::Panel object.

A standard set of options are recognized.  See L<OPTIONS>.

=back

=head2 OBJECT METHODS

Once a glyph is created, it responds to a large number of methods.  In
this section, these methods are grouped into related categories.

Retrieving glyph context:

=over 4

=item $factory = $glyph-E<gt>factory

Get the Bio::Graphics::Glyph::Factory associated with this object.
This cannot be changed once it is set.

=item $panel = $glyph-E<gt>panel

Get the Bio::Graphics::Panel associated with this object.  This cannot
be changed once it is set.

=item $feature = $glyph-E<gt>feature

Get the sequence feature associated with this object.  This cannot be
changed once it is set.

=item $feature = $glyph-E<gt>add_feature(@features)

Add the list of features to the glyph, creating subparts.  This is
most common done with the track glyph returned by
Ace::Graphics::Panel-E<gt>add_track().

=item $feature = $glyph-E<gt>add_group(@features)

This is similar to add_feature(), but the list of features is treated
as a group and can be configured as a set.

=back

Retrieving glyph options:

=over 4

=item $fgcolor = $glyph-E<gt>fgcolor

=item $bgcolor = $glyph-E<gt>bgcolor

=item $fontcolor = $glyph-E<gt>fontcolor

=item $fontcolor = $glyph-E<gt>font2color

=item $fillcolor = $glyph-E<gt>fillcolor

These methods return the configured foreground, background, font,
alternative font, and fill colors for the glyph in the form of a
GD::Image color index.

=item $color = $glyph-E<gt>tkcolor

This method returns a color to be used to flood-fill the entire glyph
before drawing (currently used by the "track" glyph).

=item $width = $glyph-E<gt>width([$newwidth])

Return the width of the glyph, not including left or right padding.
This is ordinarily set internally based on the size of the feature and
the scale of the panel.

=item $width = $glyph-E<gt>layout_width

Returns the width of the glyph including left and right padding.

=item $width = $glyph-E<gt>height

Returns the height of the glyph, not including the top or bottom
padding.  This is calculated from the "height" option and cannot be
changed.


=item $font = $glyph-E<gt>font

Return the font for the glyph.

=item $option = $glyph-E<gt>option($option)

Return the value of the indicated option.

=item $index = $glyph-E<gt>color($color)

Given a symbolic or #RRGGBB-form color name, returns its GD index.

=item $level = $glyph-E<gt>level

The "level" is the nesting level of the glyph.
Groups are level -1, top level glyphs are level 0,
subparts (e.g. exons) are level 1 and so forth.

=back

Setting an option:

=over 4

=item $glyph-E<gt>configure(-name=E<gt>$value)

You may change a glyph option after it is created using set_option().
This is most commonly used to configure track glyphs.

=back

Retrieving information about the sequence:

=over 4

=item $start = $glyph-E<gt>start

=item $end   = $glyph-E<gt>end

These methods return the start and end of the glyph in base pair
units.

=item $offset = $glyph-E<gt>offset

Returns the offset of the segment (the base pair at the far left of
the image).

=item $length = $glyph-E<gt>length

Returns the length of the sequence segment.

=back


Retrieving formatting information:

=over 4

=item $top = $glyph-E<gt>top

=item $left = $glyph-E<gt>left

=item $bottom = $glyph-E<gt>bottom

=item $right = $glyph-E<gt>right

These methods return the top, left, bottom and right of the glyph in
pixel coordinates.

=item $height = $glyph-E<gt>height

Returns the height of the glyph.  This may be somewhat larger or
smaller than the height suggested by the GlyphFactory, depending on
the type of the glyph.

=item $scale = $glyph-E<gt>scale

Get the scale for the glyph in pixels/bp.

=item $height = $glyph-E<gt>labelheight

Return the height of the label, if any.

=item $label = $glyph-E<gt>label

Return a human-readable label for the glyph.

=back

These methods are called by Bio::Graphics::Track during the layout
process:

=over 4

=item $glyph-E<gt>move($dx,$dy)

Move the glyph in pixel coordinates by the indicated delta-x and
delta-y values.

=item ($x1,$y1,$x2,$y2) = $glyph-E<gt>box

Return the current position of the glyph.

=back

These methods are intended to be overridden in subclasses:

=over 4

=item $glyph-E<gt>calculate_height

Calculate the height of the glyph.

=item $glyph-E<gt>calculate_left

Calculate the left side of the glyph.

=item $glyph-E<gt>calculate_right

Calculate the right side of the glyph.

=item $glyph-E<gt>draw($gd,$left,$top)

Optionally offset the glyph by the indicated amount and draw it onto
the GD::Image object.


=item $glyph-E<gt>draw_label($gd,$left,$top)

Draw the label for the glyph onto the provided GD::Image object,
optionally offsetting by the amounts indicated in $left and $right.

=back

These methods are useful utility routines:

=over 4

=item $pixels = $glyph-E<gt>map_pt($bases);

Map the indicated base position, given in base pair units, into
pixels, using the current scale and glyph position.

=item $glyph-E<gt>filled_box($gd,$x1,$y1,$x2,$y2)

Draw a filled rectangle with the appropriate foreground and fill
colors, and pen width onto the GD::Image object given by $gd, using
the provided rectangle coordinates.

=item $glyph-E<gt>filled_oval($gd,$x1,$y1,$x2,$y2)

As above, but draws an oval inscribed on the rectangle.

=back

=head2 OPTIONS

The following options are standard among all Glyphs.  See individual
glyph pages for more options.

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

  -strand_arrow Whether to indicate            undef (false)
                 strandedness

  -label        Whether to draw a label	       undef (false)

  -description  Whether to draw a description  undef (false)

  -sort_order   Specify layout sort order      "default"

  -always_sort  Sort even when bumping is off  undef (false)

  -bump_limit   Maximum number of levels to bump undef (unlimited)

For glyphs that consist of multiple segments, the B<-connector> option
controls what's drawn between the segments.  The default is undef (no
connector).  Options include:

   "hat"     an upward-angling conector
   "solid"   a straight horizontal connector
   "quill"   a decorated line with small arrows indicating strandedness
             (like the UCSC Genome Browser uses)
   "dashed"  a horizontal dashed line.  

The B<-connector_color> option controls the color of the connector, if
any.

The label is printed above the glyph.  You may pass an anonymous
subroutine to B<-label>, in which case the subroutine will be invoked
with the feature as its single argument.  and is expected to return
the string to use as the description.  If you provide the numeric
value "1" to B<-description>, the description will be read off the
feature's seqname(), info() and primary_tag() methods will be called
until a suitable name is found.  To create a label with the
text "1", pass the string "1 ".  (A 1 followed by a space).

The description is printed below the glyph.  You may pass an anonymous
subroutine to B<-description>, in which case the subroutine will be
invoked with the feature as its single argument and is expected to
return the string to use as the description.  If you provide the
numeric value "1" to B<-description>, the description will be read off
the feature's source_tag() method.  To create a description with the
text "1", pass the string "1 ".  (A 1 followed by a space).

In the case of ACEDB Ace::Sequence feature objects, the feature's
info(), Brief_identification() and Locus() methods will be called to
create a suitable description.

The B<-strand_arrow> option, if true, requests that the glyph indicate
which strand it is on, usually by drawing an arrowhead.  Not all
glyphs will respond to this request.  For historical reasons,
B<-stranded> is a synonym for this option.

By default, features are drawn with a layout based only on the
position of the feature, assuring a maximal "packing" of the glyphs
when bumped.  In some cases, however, it makes sense to display the
glyphs sorted by score or some other comparison, e.g. such that more
"important" features are nearer the top of the display, stacked above
less important features.  The -sort_order option allows a few
different built-in values for changing the default sort order (which
is by "left" position): "low_score" (or "high_score") will cause
features to be sorted from lowest to highest score (or vice versa).
"left" (or "default") and "right" values will cause features to be
sorted by their position in the sequence.  "longer" (or "shorter")
will cause the longest (or shortest) features to be sorted first, and
"strand" will cause the features to be sorted by strand: "+1"
(forward) then "0" (unknown, or NA) then "-1" (reverse).  Lastly,
"name" will sort features alphabetically by their display_name()
attribute.

In all cases, the "left" position will be used to break any ties.  To
break ties using another field, options may be strung together using a
"|" character; e.g. "strand|low_score|right" would cause the features
to be sorted first by strand, then score (lowest to highest), then by
"right" position in the sequence.  Finally, a subroutine coderef can
be provided, which should expect to receive two feature objects (via
the special sort variables $a and $b), and should return -1, 0 or 1
(see Perl's sort() function for more information); this subroutine
will be used without further modification for sorting.  For example,
to sort a set of database search hits by bits (stored in the features'
"score" fields), scaled by the log of the alignment length (with
"left" position breaking any ties):

  sort_order = sub { ( $b->score/log($b->length)
                                      <=>
                       $a->score/log($a->length) )
                                      ||
                     ( $a->start <=> $b->start )
                   }

The -always_sort option, if true, will sort features even if bumping
is turned off.  This is useful if you would like overlapping features
to stack in a particular order.  Features towards the end of the list
will overlay those towards the beginning of the sort order.

=head1 SUBCLASSING Bio::Graphics::Glyph

By convention, subclasses are all lower-case.  Begin each subclass
with a preamble like this one:

 package Bio::Graphics::Glyph::crossbox;

 use strict;
 use vars '@ISA';
 @ISA = 'Bio::Graphics::Glyph';

Then override the methods you need to.  Typically, just the draw()
method will need to be overridden.  However, if you need additional
room in the glyph, you may override calculate_height(),
calculate_left() and calculate_right().  Do not directly override
height(), left() and right(), as their purpose is to cache the values
returned by their calculating cousins in order to avoid time-consuming
recalculation.

A simple draw() method looks like this:

 sub draw {
  my $self = shift;
  $self->SUPER::draw(@_);
  my $gd = shift;

  # and draw a cross through the box
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  my $fg = $self->fgcolor;
  $gd->line($x1,$y1,$x2,$y2,$fg);
  $gd->line($x1,$y2,$x2,$y1,$fg);
 }

This subclass draws a simple box with two lines criss-crossed through
it.  We first call our inherited draw() method to generate the filled
box and label.  We then call calculate_boundaries() to return the
coordinates of the glyph, disregarding any extra space taken by
labels.  We call fgcolor() to return the desired foreground color, and
then call $gd-E<gt>line() twice to generate the criss-cross.

For more complex draw() methods, see Bio::Graphics::Glyph::transcript
and Bio::Graphics::Glyph::segments.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::DB::GFF::Feature>,
L<Ace::Sequence>,
L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>,
L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::dna>,
L<Bio::Graphics::Glyph::graded_segments>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::wormbase_transcript>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
