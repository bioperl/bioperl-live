package Bio::Graphics::Glyph::arrow;
# package to use for drawing an arrow

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::generic;
@ISA = 'Bio::Graphics::Glyph::generic';

my %UNITS = (n => 1e-12,
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
sub draw {
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

  $trunc_left = 0  if $self->no_trunc;
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

    my $relative = $self->option('relative_coords');
    my $relative_coords_offset = $self->option('relative_coords_offset');
    $relative_coords_offset = 1 unless ($relative_coords_offset);

    my $start    = $relative ? $relative_coords_offset : $self->feature->start-1;
    my $stop     = $start + $self->feature->length - 1;

    my $offset   = $relative ? ($self->feature->start - $relative_coords_offset) : 0;
    my $reversed = exists $self->{flip} || ($relative && $self->feature->strand < 0);

    my $unit_label   = $self->option('units') || '';
    my $unit_divider = $self->option('unit_divider') || 1;

    my $units      = $self->calculate_units($start/$unit_divider,$self->feature->length/$unit_divider);
    my $divisor    = $UNITS{$units} || 1;

    $divisor *= $unit_divider;

    my $format     = min($self->feature->length,$self->panel->length)/$divisor > 10
      ? "%d$units%s" : "%.6g$units%s";

    my $scale  = $self->option('scale') || 1;  ## Does the user want to override the internal scale?

    my $model  = sprintf("$format ",$stop/($divisor*$scale),$unit_label);
    my $minlen = $width * length($model);

    my ($major_interval,$minor_interval) = $self->panel->ticks(($stop-$start+1)/$unit_divider,$minlen);

    my $left  = $sw ? $x1+$height : $x1;
    my $right = $ne ? $x2-$height : $x2;

    # adjust for portions of arrow that are outside panel
    $start += $self->panel->start - $self->feature->start
      if $self->feature->start < $self->panel->start;
    $stop  -= $self->feature->end - $self->panel->end
      if $self->feature->end   > $self->panel->end;

    my $first_tick = $major_interval * int(0.5 + $start/$major_interval);
    my $last_tick  = $major_interval * int(0.5 + $stop/$major_interval);

    for (my $i = $first_tick; $i <= $last_tick; $i += $major_interval) {

      my $tickpos = $dx + ($reversed ? $self->map_pt($stop - $i + $offset)
	                             : $self->map_pt($i + $offset));
      next if $tickpos < $left or $tickpos > $right;

      $gd->line($tickpos,$center-$a2,$tickpos,$center+$a2,$fg);
      my $label = $scale ? $i / $scale : $i;
      my $scaled = $label/$divisor;
      $label = sprintf($format,$scaled,$unit_label);

      my $middle = $tickpos - (length($label) * $width)/2;
      next if $middle < $left or $middle > $right;

      $gd->string($font,$middle,$center+$a2-1,$label,$font_color)
        unless ($self->option('no_tick_label'));
    }

    if ($self->option('tick') >= 2) {

      $first_tick = $minor_interval * int(0.5 + $start/$minor_interval);
      $last_tick  = $minor_interval * int(0.5 + $stop/$minor_interval);

      my $a4 = $self->height/4;
      for (my $i = $first_tick; $i <= $last_tick; $i += $minor_interval) {
	my $tickpos = $dx + ($reversed ? $self->map_pt($stop - $i + $offset)
	                               : $self->map_pt($i + $offset));
	next if $tickpos < $left or $tickpos > $right;

	$gd->line($tickpos,$center-$a4,$tickpos,$center+$a4,$fg);
      }
    }
  }

  # add a label if requested
  $self->draw_label($gd,$dx,$dy)       if $self->option('label');
  $self->draw_description($gd,$dx,$dy) if $self->option('description');
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

Ace::Graphics::Glyph::arrow - The "arrow" glyph

=head1 SYNOPSIS

  See L<Ace::Graphics::Panel> and L<Ace::Graphics::Glyph>.

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

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -tick       Whether to draw major             0
              and minor ticks.
	      0 = no ticks
	      1 = major ticks
	      2 = minor ticks

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
