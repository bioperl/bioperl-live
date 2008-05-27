package Bio::Graphics::Glyph::arrow;
# package to use for drawing an arrow

# $Id$
# Non object-oriented utilities used here-and-there in Bio::Graphics modules

=head1 NAME

Bio::Graphics::Glyph::arrow - the "arrow" glyph

=cut

use strict;
use Bio::Coordinate::Pair;
use Bio::Location::Simple;
use base qw(Bio::Graphics::Glyph::generic);

my %UNITS = (p => 1e-12,
	     n => 1e-9,
	     u => 1e-6,
	     m => 0.001,
	     c => 0.01,
	     k => 1000,
	     M => 1_000_000,
	     G => 1_000_000_000);

sub pad_bottom {
  my $self = shift;
  my $val = $self->SUPER::pad_bottom(@_);
  $val += $self->font->height if $self->option('tick');
  $val;
}

# override draw method
sub draw_component {
  my $self = shift;
  my $parallel = $self->option('parallel');
  $parallel = 1 unless defined $parallel;
  $self->draw_parallel(@_) if $parallel;
  $self->draw_perpendicular(@_) unless $parallel;
}

sub draw_perpendicular {
  my $self = shift;
  my $gd = shift;
  my ($dx,$dy) = @_;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $ne = $self->option('northeast');
  my $sw = $self->option('southwest');
  $ne = $sw = 1 unless defined($ne) || defined($sw);

  # draw a perpendicular arrow at position indicated by $x1
  my $fg = $self->set_pen;
  my $a2 = ($y2-$y1)/4;

  my @positions = $x1 == $x2 ? ($x1) : ($x1,$x2);
  for my $x (@positions) {
    if ($ne) {
      $gd->line($x,$y1,$x,$y2,$fg);
      $gd->line($x-$a2,$y1+$a2,$x,$y1,$fg);
      $gd->line($x+$a2,$y1+$a2,$x,$y1,$fg);
    }
    if ($sw) {
      $gd->line($x,$y1,$x,$y2,$fg);
      $gd->line($x-$a2,$y2-$a2,$x,$y2,$fg);
      $gd->line($x+$a2,$y2-$a2,$x,$y2,$fg);
    }
  }

  # add a label if requested
  $self->draw_label($gd,$dx,$dy) if $self->option('label');  # this draws the label aligned to the left
}

sub draw_parallel {
  my $self = shift;
  my $gd = shift;
  my ($dx,$dy) = @_;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $fg = $self->set_pen;
  my $a2 = ($self->height)/2;
  my $center = $y1+$a2;

  my $trunc_left  = $x1 < $self->panel->left;
  my $trunc_right = $x2 > $self->panel->right;
  $x1 = $self->panel->left  if $trunc_left;
  $x2 = $self->panel->right if $trunc_right;

#  warn $self->feature,": x1=$x1, x2=$x2, start=$self->{start},end=$self->{end}, strand=$self->{strand}";
#  warn join ' ',%$self;

  $trunc_left  = 0 if $self->no_trunc;
  $trunc_right = 0 if $self->no_trunc;

  my ($sw,$ne,$base_w,$base_e) = $self->arrowheads;
  $gd->line($x1,$center,$x2,$center,$fg);
  $self->arrowhead($gd,$x1,$center,$a2,-1) if $sw && !$trunc_left;  # west arrow
  $self->arrowhead($gd,$x2,$center,$a2,+1) if $ne && !$trunc_right; # east arrow
  $gd->line($x1,$center-$a2,$x1,$center+$a2,$fg) if $base_w && !$trunc_left;  #west base
  $gd->line($x2,$center-$a2,$x2,$center+$a2,$fg) if $base_e && !$trunc_right; #east base

  # turn on ticks
  if ($self->option('tick')) {
    local $^W = 0;  # dumb uninitialized variable warning
    my $font       = $self->font;
    my $width      = $font->width;
    my $font_color = $self->fontcolor;
    my $height     = $self->height;

    my $relative   = $self->option('relative_coords');
    my $flipped    = $self->{flip};
    my $end        = $self->panel->end + 1;

    my $tickwidth  = $self->option('tickwidth'); $tickwidth = $self->linewidth unless defined $tickwidth;
    my $tickcolor  = $self->color($self->option('tickcolor') || $self->option('fgcolor'));
    my $tickpen    = $self->set_pen($tickwidth, $tickcolor);

    my $relative_coords_offset = $self->option('relative_coords_offset');
    $relative_coords_offset    = 1 unless defined $relative_coords_offset;

    my $start    = $relative ? $relative_coords_offset : $self->feature->start-1;
    my $stop     = $start + $self->feature->length - 1;

    my $map = Bio::Coordinate::Pair->new(-in  => Bio::Location::Simple->new( -seq_id => "rel",
									     -start => $start,
									     -end => $stop,
									     -strand => 1,
									     ),
					 -out => Bio::Location::Simple->new( -seq_id => "abs",
									     -start => $self->feature->start,
									     -end => $self->feature->end,
									     -strand => $self->feature->strand,
									     ),
					 ) if $relative;

    my $unit_label     = $self->option('units')        || '';
    my $unit_divider   = $self->option('unit_divider') || 1;
    my $units_in_label = $self->option('units_in_label');

    my $units      = $self->calculate_units($start/$unit_divider,$self->feature->length/$unit_divider);
    my $divisor    = $UNITS{$units} || 1;

    $divisor *= $unit_divider;

    my $format     = min($self->feature->length,$self->panel->length)/$divisor > 10
      ? "%d" : "%.6g";

    $format .= "$units%s" unless $units_in_label;

    my $scale  = $self->option('scale') || 1;  ## Does the user want to override the internal scale?

    my $model  = sprintf("$format ",$stop/($divisor*$scale),$unit_label);
    $model     = "-$model" if $start < 0;

    my $minlen = $width * length($model);# * 1.5;

    my ($major_interval,$minor_interval) = $self->panel->ticks(($stop-$start+1)/$unit_divider,$minlen);

    my $left  = $sw ? $x1+$height : $x1;
    my $right = $ne ? $x2-$height : $x2;

    # adjust for portions of arrow that are outside panel
    if ($relative && $self->feature->strand == -1) {
	$start += $self->feature->end - $self->panel->end if $self->feature->end > $self->panel->end;
	$stop -= $self->panel->start - $self->feature->start if $self->feature->start < $self->panel->start;
    } else {
	$start += $self->panel->start - $self->feature->start
	    if $self->feature->start < $self->panel->start;
	$stop  -= $self->feature->end - $self->panel->end
	    if $self->feature->end   > $self->panel->end;
    }
	
    my $first_tick = $major_interval * int($start/$major_interval);
    my $last_tick  = $major_interval * int(($stop+2)/$major_interval);

    my $label_intervals = $self->label_intervals;
    my $interval_width  = $major_interval * $self->scale/2;
    my %drewit;

    for (my $i = $first_tick; $i <= $last_tick; $i += $major_interval) {
      my $abs = $i;
      if ($relative) {
	  $abs = $map->map( Bio::Location::Simple->new(-seq_id => "rel",
						       -start  => $i,
						       -end   => $i,
						       -strand => 1,
						       )
			    )->match;
	  next unless $abs;
	  $abs = $abs->start;
      }

      $abs = $end - $abs + 1 if $flipped;

      my $tickpos = int $dx + $self->map_pt($abs);
      next if $tickpos < $x1 || $tickpos > $x2;
      $drewit{$tickpos}++;

      $gd->line($tickpos,$center-$a2,$tickpos,$center+$a2,$tickpen)
	unless $tickpos < $left or $tickpos > $right;

      my $label = $scale ? $i / $scale : $i;
      my $scaled = $label/$divisor;
      $label = sprintf($format,$scaled,$unit_label);

      my $label_len = length($label) * $width;

      my $middle = $tickpos - $label_len/2;
      $middle   += $interval_width if $label_intervals;

      $gd->string($font,$middle,$center+$a2-1,$label,$font_color)
        unless ($self->option('no_tick_label') || $middle > $x2);
    }

    if ($self->option('tick') >= 2) {

      $first_tick = $minor_interval * int($start/$minor_interval);
      $last_tick  = $minor_interval * int(($stop+2)/$minor_interval);

      my $a4 = $self->height/4;
      for (my $i = $first_tick; $i <= $last_tick; $i += $minor_interval) {
	  my $abs = $i;
	  if ($relative) {
	      $abs = $map->map( Bio::Location::Simple->new(-seq_id => "rel",
							   -start  => $i,
							   -end    => $i,
							   -strand => 1,
							   )
				)->match;
	      next unless $abs;
	      $abs = $abs->start;
	  }
	  $abs = $end - $abs if $flipped;

	  my $tickpos = int $dx + $self->map_pt($abs);
	  next if $tickpos < $left-1 or $tickpos > $right+1;
	  next if $drewit{$tickpos} || $drewit{$tickpos-1} || $drewit{$tickpos+1}; # prevent roundoff errors from appearing

	  $gd->line($tickpos,$center-$a4,$tickpos,$center+$a4,$tickpen);
      }
    }
  }

  # add a label if requested
  $self->draw_label($gd,$dx,$dy)       if $self->option('label');
  $self->draw_description($gd,$dx,$dy) if $self->option('description');
}

sub label {
  my $self  = shift;
  my $label = $self->SUPER::label(@_);
  return $label unless $self->option('units_in_label');
  my $unit_divider = $self->option('unit_divider') || 1;
  my $unit_label   = $self->option('units')        || '';
  my $start        = $self->feature->start-1;
  my $units        = $self->calculate_units($start/$unit_divider,$self->feature->length/$unit_divider);
  return $label . " ($units$unit_label)";
}

sub label_intervals {
  return shift->option('label_intervals');
}

sub arrowheads {
  my $self = shift;
  my ($ne,$sw,$base_e,$base_w);
  if ($self->option('double')) {
    $ne = $sw = 1;
  } else {
    $ne   = $self->option('northeast') || $self->option('east');
    $sw   = $self->option('southwest') || $self->option('west');
  }
  # otherwise use strandedness to define the arrow
  unless (defined($ne) || defined($sw)) {
    # turn on both if neither specified
    $ne = 1 if $self->feature->strand > 0;
    $sw = 1 if $self->feature->strand < 0;
    ($ne,$sw) = ($sw,$ne) if $self->{flip};
  }
  return ($sw,$ne,0,0) unless $self->option('base');
  return ($sw,$ne,
	  (!$sw && $self->feature->start>= $self->panel->start),
	  (!$ne && $self->feature->end  <= $self->panel->end));
}

sub no_trunc { 0; }

sub calculate_units {
  my $self   = shift;
  my ($start,$length) = @_;
  return 'G' if $length >= 1e9;
  return 'M' if $length >= 1e6;
  return 'k' if $length >= 1e3;
  return ''  if $length >= 1;
  return 'c' if $length >= 1e-2;
  return 'm' if $length >= 1e-3;
  return 'u' if $length >= 1e-6;
  return 'n' if $length >= 1e-9;
  return 'p';
}

sub min { $_[0]<$_[1] ? $_[0] : $_[1] }

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::arrow - The "arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws arrows.  Depending on options, the arrows can be
labeled, be oriented vertically or horizontally, or can contain major
and minor ticks suitable for use as a scale.

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

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -tick       Whether to draw major             0
              and minor ticks.
	      0 = no ticks
	      1 = major ticks
	      2 = minor ticks

  -tickcolor  Color to use for tick marks       fgcolor

  -tickwidth  Line width to use for ticks       linewidth

  -parallel   Whether to draw the arrow         1 (true)
	      parallel to the sequence
	      or perpendicular to it.

  -northeast  Force a north or east             1 (true)
	      arrowhead(depending 
	      on orientation)

  -east       synonym of above

  -southwest  Force a south or west             1 (true)
	      arrowhead(depending 
	      on orientation)

  -west       synonym of above

  -double     force-doubleheaded arrow          0 (false)

  -base       Draw a vertical base at the       0 (false)
              non-arrowhead side

  -scale      Reset the labels on the arrow     0 (false)
              to reflect an externally 
              established scale.

  -arrowstyle "regular" to create a simple      regular
              arrowhead.  "filled" to create
              a thick filled arrowhead

  -relative_coords 
                 use relative coordinates       0 (false)
                 for scale

  -relative_coords_offset 
                 set the relative offset        1 
                 for scale

  -label_intervals                              0 (false)
              Put the numeric labels on the
              intervals between the ticks 
              rather than on the ticks
              themselves.

  -units      add units to the tick labels      none
              e.g. bp

  -unit_divider                                 1
              divide tick labels by the
              indicated amount prior to
              displaying (use, for example
              if you want to display in
              cR units)

Set -parallel to 0 (false) to display a point-like feature such as a
polymorphism, or to indicate an important location.  If the feature
start == end, then the glyph will draw a single arrow at the
designated location:

       ^
       |

Otherwise, there will be two arrows at the start and end:

       ^              ^
       |              |

Scale: Pass in a externally established scale to reset the labels on
the arrow.  This is particularly useful for manually constructed
images where the founding parameters of the panel are not 1-based.
For example, a genetic map interval ranging from 0.1 - 0.3 can be
constructed by first multiplying every value by 100. Passing

  arrow(-scale=>100);

will draw tick marks labelled appropriately to your external scale.

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

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
