package Bio::Graphics::Glyph;
use GD;

use strict;
use Carp 'croak';
use constant BUMP_SPACING => 2; # vertical distance between bumped glyphs
use vars '$VERSION';
$VERSION = '1.02';

my %LAYOUT_COUNT;

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

  my $self = bless {},$class;
  $self->{feature} = $feature;
  $self->{factory} = $factory;
  $self->{level}   = $level;
  $self->{top} = 0;

  my @subglyphs;
  my @subfeatures = $self->subseq($feature);

  if (@subfeatures) {

    # dynamic glyph resolution
    @subglyphs = sort { $a->left  <=> $b->left }  $factory->make_glyph($level+1,@subfeatures);

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
      my $right        = (sort { $b<=>$a } map {$_->right} @subglyphs)[0];
      my $w            = $right - $self->{left} + 1;
      $self->{width}   = $w if !defined($self->{width}) || $w > $self->{width};
  }

  #Handle glyphs that don't actually fill their space, but merely mark a point.
  #They need to have their collision bounds altered.  We will (for now)
  #hard code them to be in the center of their feature.
  $self->{point} = $arg{-point} ? $self->height : undef;
  if($self->option('point')){
    my ($left,$right) = $factory->map_pt($self->start,$self->stop);
    my $center = int(($left+$right)/2);

    $self->{width} = $self->height;
    $self->{left}  = $center - ($self->{width});
    $self->{right} = $center + ($self->{width});
  }

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
  $self->{start} = $self->{feature}->start;

  # handle the case of features whose endpoints are undef
  # (this happens with wormbase clones where one or more clone end is not defined)
  # in this case, we set the start to one minus the beginning of the panel
  $self->{start} = $self->panel->offset - 1 unless defined $self->{start};

  return $self->{start};
}
sub stop    {
  my $self = shift;
  return $self->{stop} if exists $self->{stop};
  $self->{stop} = $self->{feature}->end;

  # handle the case of features whose endpoints are undef
  # (this happens with wormbase clones where one or more clone end is not defined)
  # in this case, we set the start to one plus the end of the panel
  $self->{stop} = $self->panel->offset + $self->panel->length + 1 unless defined $self->{stop};

  return $self->{stop}
}
sub end     { shift->stop }
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
#  return $self->{cache_left} if exists $self->{cache_left};
#  $self->{cache_left} = $self->{left} - $self->pad_left;
  return $self->{left} - $self->pad_left;
}
sub right {
  my $self = shift;
#  return $self->{cache_right} if exists $self->{cache_right};
#  $self->{cache_right} = $self->left + $self->layout_width - 1;
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
#  return $self->{layout_width} ||= $self->width + $self->pad_left + $self->pad_right;
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
  my ($x1,$y1,$x2,$y2) = @_;

  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;
  my $linewidth = $self->option('linewidth') || 1;

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


# return boxes surrounding each part
sub boxes {
  my $self = shift;
  my ($left,$top) = @_;
  $top  += 0; $left += 0;
  my @result;

  $self->layout;
  for my $part ($self->parts) {
    if (eval{$part->feature->primary_tag} eq 'group') {
      push @result,$part->boxes($left+$self->left,$top+$self->top);
    } else {
      my ($x1,$y1,$x2,$y2) = $part->box;
      push @result,[$part->feature,$x1,$top+$self->top+$y1,$x2,$top+$self->top+$y2];
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
#              +2   simple (hyper) bump up
#              -2   simple (hyper) bump down
sub bump {
  my $self = shift;
  return $self->option('bump');
}

# we also look for the "color" option for Ace::Graphics compatibility
sub fgcolor {
  my $self = shift;
  my $index = $self->option('fgcolor') || $self->option('color') || return 0;
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
  my $index = $self->option('bgcolor') || $self->option('fillcolor') || return 0;
  $self->factory->translate_color($index);
}
sub font {
  shift->option('font');
}
sub fontcolor {
  my $self = shift;
  $self->color('fontcolor') || $self->fgcolor;
}
sub font2color {
  my $self = shift;
  $self->color('font2color') || $self->fontcolor;
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

# handle collision detection
sub layout {
  my $self = shift;
  return $self->{layout_height} if exists $self->{layout_height};

  (my @parts = $self->parts)
    || return $self->{layout_height} = $self->height + $self->pad_top + $self->pad_bottom;

  my $bump_direction = $self->bump;

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

  if (abs($bump_direction) <= 1) {  # original bump algorithm

    my %occupied; # format of occupied: key={top,bottom}, value=right
    for my $g (sort { $a->left <=> $b->left } @parts) {

      my $pos = 0;
      my $left   = $g->left;
      my $right  = $g->right;
      my $height = $g->{layout_height};

      while (1) {
	# look for collisions
	my $bottom = $pos + $height;

	my $collision;
	for my $key (keys %occupied) {
	  my ($oldtop,$oldbottom) = split /,/,$key;
	  my $oldright = $occupied{$key};
	  next if $oldright+2  < $left;
	  next if $oldbottom   < $pos;
	  next if $oldtop      > $bottom;
	  $collision = [$oldtop,$oldbottom,$oldright];
	  last;
	}
	last unless $collision;

	if ($bump_direction > 0) {
	  $pos += $collision->[1]-$collision->[0] + BUMP_SPACING;    # collision, so bump

	} else {
	  $pos -= BUMP_SPACING;
	}

      }

      $g->move(0,$pos);
      my $key = join ',',$g->top,$g->bottom;
      $occupied{$key} = $right if !exists $occupied{$key} or $occupied{$key} < $right;
    }
  }

  else {  # abs(bump) >= 2 -- simple bump algorithm
    my $pos = 0;
    my $last;
    for my $g (sort { $a->left <=> $b->left } @parts) {
      next if !defined($last);
      $pos += $bump_direction > 0 ? $last->{layout_height} + BUMP_SPACING 
                                  : - ($g->{layout_height}+BUMP_SPACING);
      $g->move(0,$pos);
    } continue {
      $last = $g;
    }
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

sub draw {
  my $self = shift;
  my $gd = shift;
  my ($left,$top,$partno,$total_parts) = @_;

  local($self->{partno},$self->{total_parts});
  @{$self}{qw(partno total_parts)} = ($partno,$total_parts);

  my $connector =  $self->connector;
  if (my @parts = $self->parts) {
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
    $self->_connector($gd,$dx,$dy,$x1,$xt,$x1,$xb,$xl,$xt,$xr,$xb);
    my ($xl2,$xt2,$xr2,$xb2) = $parts[-1]->bounds;
    $self->_connector($gd,$dx,$dy,$parts[-1]->bounds,$x2,$xt2,$x2,$xb2) if $xr2 >= $xr;
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
  #    return unless $right-$left > 1;

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

  if ($self->option('strand_arrow')) {
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
  $factory->set_option(label => 1);
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

  -connector    Connector type                 0 (false)

  -connector_color
                Connector color                black

  -strand_arrow Whether to indicate            0 (false)
                 strandedness

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

For glyphs that consist of multiple segments, the -connector option
controls what's drawn between the segments.  The default is 0 (no
connector).  Options include "hat", an upward-angling conector,
"solid", a straight horizontal connector, and "dashed", for a
horizontal dashed line.  The -connector_color option controls the
color of the connector, if any.

The label is printed above the glyph.  You may pass an anonymous
subroutine to -label, in which case the subroutine will be invoked
with the feature as its single argument.  The subroutine must return a
string to render as the label.  Otherwise, you may return the number
"1", in which case the feature's info(), seqname() and primary_tag()
methods will be called (in that order) until a suitable name is found.

The description is printed below the glyph.  You may pass an anonymous
subroutine to -label, in which case the subroutine will be invoked
with the feature as its single argument.  The subroutine must return a
string to render as the label.  Otherwise, you may return the number
"1", in which case the feature's source_tag() method will be invoked.

In the case of ACEDB Ace::Sequence feature objects, the feature's
info(), Brief_identification() and Locus() methods will be called to
create a suitable description.

The -strand_arrow option, if true, requests that the glyph indicate
which strand it is on, usually by drawing an arrowhead.  Not all
glyphs can respond appropriately to this request.

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
