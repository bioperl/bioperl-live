package Bio::Graphics::Glyph::xyplot;

use strict;
use Bio::Graphics::Glyph::segments;
use vars '@ISA';
use GD 'gdTinyFont';

@ISA = 'Bio::Graphics::Glyph::segments';

use constant DEFAULT_POINT_RADIUS=>1;

my %SYMBOLS = (
	       triangle => \&draw_triangle,
	       square   => \&draw_square,
	       disc     => \&draw_disc,
	       point    => \&draw_point,
	      );

sub point_radius {
  shift->option('point_radius') || DEFAULT_POINT_RADIUS;
}

sub pad_top { 0 }

sub draw {
  my $self = shift;
  my ($gd,$dx,$dy) = @_;
  my ($left,$top,$right,$bottom) = $self->calculate_boundaries($dx,$dy);

  my @parts = $self->parts;
  return $self->SUPER::draw(@_) unless @parts > 0;

  # figure out the scale and such like
  my $max_score = $self->option('max_score');
  my $min_score = $self->option('min_score');

  unless (defined $max_score && defined $min_score) {
    my $first = $parts[0];
    $max_score = $min_score = eval { $first->feature->score} || 0;
    for my $part (@parts) {
      my $s = eval { $part->feature->score };
      next unless defined $s;
      $max_score = $s if $s > $max_score;
      $min_score = $s if $s < $min_score;
    }
  }

  # if a scale is called for, then we adjust the max and min to be even
  # multiples of a power of 10.
  if ($self->option('scale')) {
    $max_score = max10($max_score);
    $min_score = min10($min_score);
  }

  my $height = $self->option('height');
  my $scale  = $max_score > $min_score ? $height/($max_score-$min_score)
                                       : 1;
  my $x = $dx;
  my $y = $dy + $self->top + $self->pad_top;

  # now seed all the parts with the information they need to draw their positions
  foreach (@parts) {
    my $s = eval {$_->feature->score};
    next unless defined $s;
    my $position      = ($s-$min_score) * $scale;
    $_->{_y_position} = $bottom - $position;
  }

  my $type = $self->option('graph_type');
  $self->_draw_histogram($gd,$x,$y)  if $type eq 'histogram';
  $self->_draw_boxes($gd,$x,$y)      if $type eq 'boxes';
  $self->_draw_line ($gd,$x,$y)      if $type eq 'line'
                                       or $type eq 'linepoints';
  $self->_draw_points($gd,$x,$y)     if $type eq 'points'
                                       or $type eq 'linepoints';

  $self->_draw_scale($gd,$scale,$min_score,$max_score,$dx,$dy)      if $self->option('scale');
}

sub log10 { log(shift)/log(10) }
sub max10 {
  my $a = shift;
  $a = 1 if $a <= 0;
  my $l=int(log10($a)); 
  $l = 10**$l; 
  my $r = $a/$l; 
  return $r*$l if int($r) == $r;
  return $l*int(($a+$l)/$l);
}
sub min10 {
  my $a = shift;
  $a = 1 if $a <= 0;
  my $l=int(log10($a));
  $l = 10**$l; 
  my $r = $a/$l; 
  return $r*$l if int($r) == $r;
  return $l*int($a/$l);
}

sub _draw_histogram {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my @parts  = $self->parts;
  my $fgcolor = $self->fgcolor;

  # draw each of the component lines of the histogram surface
  for (my $i = 0; $i < @parts; $i++) {
    my $part = $parts[$i];
    my $next = $parts[$i+1];
    my ($x1,$y1,$x2,$y2) = $part->calculate_boundaries($left,$top);
    $gd->line($x1,$part->{_y_position},$x2,$part->{_y_position},$fgcolor);
    next unless $next;
    my ($x3,$y3,$x4,$y4) = $next->calculate_boundaries($left,$top);
    if ($x2 == $x3) {# connect vertically to next level
      $gd->line($x2,$part->{_y_position},$x2,$next->{_y_position},$fgcolor); 
    } else {
      $gd->line($x2,$part->{_y_position},$x2,$y2,$fgcolor); # to bottom
      $gd->line($x2,$y2,$x3,$y2,$fgcolor);                        # to right
      $gd->line($x3,$y4,$x3,$next->{_y_position},$fgcolor);   # up
    }
  }

  # end points: from bottom to first
  my ($x1,$y1,$x2,$y2) = $parts[0]->calculate_boundaries($left,$top);
  $gd->line($x1,$y2,$x1,$parts[0]->{_y_position},$fgcolor);
  # from last to bottom
  my ($x3,$y3,$x4,$y4) = $parts[-1]->calculate_boundaries($left,$top);
  $gd->line($x4,$parts[-1]->{_y_position},$x4,$y4,$fgcolor);

  # from left to right  -- don't like this
  # $gd->line($x1,$y2,$x4,$y4,$fgcolor);

  # That's it.  Not too hard.
}

sub _draw_boxes {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my @parts  = $self->parts;
  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;
  my $height  = $self->height;

  # draw each of the component lines of the histogram surface
  for (my $i = 0; $i < @parts; $i++) {
    my $part = $parts[$i];
    my $next = $parts[$i+1];
    my ($x1,$y1,$x2,$y2) = $part->calculate_boundaries($left,$top);
    $self->filled_box($gd,$x1,$part->{_y_position},$x2,$y2,$bgcolor,$fgcolor);
    next unless $next;
    my ($x3,$y3,$x4,$y4) = $next->calculate_boundaries($left,$top);
    $gd->line($x2,$y2,$x3,$y4,$fgcolor) if $x2 < $x3;
  }

  # That's it.
}

sub _draw_line {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my @parts  = $self->parts;
  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;

  # connect to center positions of each interval
  my $first_part = shift @parts;
  my ($x1,$y1,$x2,$y2) = $first_part->calculate_boundaries($left,$top);
  my $current_x = ($x1+$x2)/2;
  my $current_y = $first_part->{_y_position};

  for my $part (@parts) {
    my ($x1,$y1,$x2,$y2) = $part->calculate_boundaries($left,$top);
    my $next_x = ($x1+$x2)/2;
    my $next_y = $part->{_y_position};
    $gd->line($current_x,$current_y,$next_x,$next_y,$fgcolor);
    ($current_x,$current_y) = ($next_x,$next_y);
  }

}

sub _draw_points {
  my $self = shift;
  my ($gd,$left,$top) = @_;
  my $symbol_name = $self->option('point_symbol') || 'point';
  my $symbol_ref  = $SYMBOLS{$symbol_name};

  my @parts   = $self->parts;
  my $bgcolor = $self->bgcolor;
  my $pr      = $self->point_radius;

  for my $part (@parts) {
    my ($x1,$y1,$x2,$y2) = $part->calculate_boundaries($left,$top);
    my $x = ($x1+$x2)/2;
    my $y = $part->{_y_position};
    $symbol_ref->($gd,$x,$y,$pr,$bgcolor);
  }
}

sub _draw_scale {
  my $self = shift;
  my ($gd,$scale,$min,$max,$dx,$dy) = @_;
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries($dx,$dy);

  my $side = $self->option('scale');
  return if $side eq 'none';
  $side   ||= 'both';

  my $fg    = $self->fgcolor;
  my $half  = ($y1+$y2)/2;

  $gd->line($x1+1,$y1,$x1+1,$y2,$fg) if $side eq 'left'  || $side eq 'both';
  $gd->line($x2-2,$y1,$x2-2,$y2,$fg) if $side eq 'right' || $side eq 'both';

  for ([$y1,$max],[$half,int(($max-$min)/2+0.5)]) {
    $gd->line($x1,$_->[0],$x1+3,$_->[0],$fg) if $side eq 'left'  || $side eq 'both';
    $gd->line($x2-4,$_->[0],$x2,$_->[0],$fg) if $side eq 'right' || $side eq 'both';
    if ($side eq 'left' or $side eq 'both') {
      $gd->string(gdTinyFont,
		  $x1 + 5,$_->[0]-(gdTinyFont->height/3),
		  $_->[1],
		  $fg);
    }
    if ($side eq 'right' or $side eq 'both') {
      $gd->string(gdTinyFont,
		  $x2-5 - (length($_->[1])*gdTinyFont->width),$_->[0]-(gdTinyFont->height/3),
		  $_->[1],
		  $fg);
    }
  }
}

# we are unbumpable!
sub bump {
  return 0;
}

sub connector {
  my $self = shift;
  my $type = $self->option('graph_type');
  return 1 if $type eq 'line' or $type eq 'linepoints';
}

sub height {
  my $self = shift;
  return $self->option('graph_height') || $self->SUPER::height;
}

sub draw_triangle {
  my ($gd,$x,$y,$pr,$color) = @_;
  my ($vx1,$vy1) = ($x-$pr,$y+$pr);
  my ($vx2,$vy2) = ($x,  $y-$pr);
  my ($vx3,$vy3) = ($x+$pr,$y+$pr);
  $gd->line($vx1,$vy1,$vx2,$vy2,$color);
  $gd->line($vx2,$vy2,$vx3,$vy3,$color);
  $gd->line($vx3,$vy3,$vx1,$vy1,$color);
}
sub draw_square {
  my ($gd,$x,$y,$pr,$color) = @_;
  $gd->line($x-$pr,$y-$pr,$x+$pr,$y-$pr,$color);
  $gd->line($x+$pr,$y-$pr,$x+$pr,$y+$pr,$color);
  $gd->line($x+$pr,$y+$pr,$x-$pr,$y+$pr,$color);
  $gd->line($x-$pr,$y+$pr,$x-$pr,$y-$pr,$color);
}
sub draw_disc {
  my ($gd,$x,$y,$pr,$color) = @_;
  $gd->arc($x,$y,$pr,$pr,0,360,$color);
}
sub draw_point {
  my ($gd,$x,$y,$pr,$color) = @_;
  $gd->setPixel($x,$y,$color);
}

sub _subseq {
  my $class   = shift;
  my $feature = shift;
  return $feature->segments                if $feature->can('segments');
  my @split = eval { my $id   = $feature->location->seq_id;
		     my @subs = $feature->location->sub_Location;
		     grep {$id eq $_->seq_id} @subs};
  return @split if @split;
  return $feature->sub_SeqFeature          if $feature->can('sub_SeqFeature');
  return;
}

sub keyglyph {
  my $self = shift;

  my $scale = 1/$self->scale;  # base pairs/pixel

  my $feature =
    Bio::Graphics::Feature->new(
				-segments=>[ [ 0*$scale,9*$scale],
					     [ 10*$scale,19*$scale],
					     [ 20*$scale, 29*$scale]
					   ],
				-name => 'foo bar',
				-strand => '+1');
  ($feature->segments)[0]->score(10);
  ($feature->segments)[1]->score(50);
  ($feature->segments)[2]->score(25);
  my $factory = $self->factory->clone;
  $factory->set_option(label => 1);
  $factory->set_option(bump  => 0);
  $factory->set_option(connector  => 'solid');
  my $glyph = $factory->make_glyph(0,$feature);
  return $glyph;
}


1;

__END__

=head1 NAME

Bio::Graphics::Glyph::xyplot - The xyplot glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing features that have a position on the
genome and a numeric value.  It can be used to represent gene
prediction scores, motif-calling scores, percent similarity,
microarray intensities, or other features that require a line plot.

The X axis represents the position on the genome, as per all other
glyphs.  The Y axis represents the score.  Options allow you to set
the height of the glyph, the maximum and minimum scores, the color of
the line and axis, and the symbol to draw.

The plot is designed to work on a single feature group that contains
subfeatures.  It is the subfeatures that carry the score
information. The best way to arrange for this is to create an
aggregator for the feature.  We'll take as an example a histogram of
repeat density in which interval are spaced every megabase and the
score indicates the number of repeats in the interval; we'll assume
that the database has been loaded in in such a way that each interval
is a distinct feature with the method name "density" and the source
name "repeat".  Furthermore, all the repeat features are grouped
together into a single group (the name of the group is irrelevant).
If you are using Bio::DB::GFF and Bio::Graphics directly, the sequence
of events would look like this:

  my $agg = Bio::DB::GFF::Aggregator->new(-method    => 'repeat_density',
                                          -sub_parts => 'density:repeat');
  my $db  = Bio::DB::GFF->new(-dsn=>'my_database',
                              -aggregators => $agg);
  my $segment  = $db->segment('Chr1');
  my @features = $segment->features('repeat_density');

  my $panel = Bio::Graphics::Panel->new;
  $panel->add_track(\@features,
                    -glyph => 'xyplot');

If you are using Generic Genome Browser, you will add this to the
configuration file:

  aggregators = repeat_density{density:repeat}
                clone alignment etc

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

  -label        Whether to draw a label	       0 (false)

  -description  Whether to draw a description  0 (false)

In addition, the alignment glyph recognizes the following
glyph-specific options:

  Option         Description                  Default
  ------         -----------                  -------

  -max_score   Maximum value of the	      Calculated
               feature's "score" attribute

  -min_score   Minimum value of the           Calculated
               feature's "score" attribute

  -graph_type  Type of graph to generate.     Histogram
               Options are: "histogram",
               "boxes", "line", "points",
               or "linepoints".

  -point_symbol Symbol to use. Options are    none
                "triangle", "square", "disc",
                "point", and "none".

  -point_radius Radius of the symbol, in      1
                pixels

  -scale        Position where the Y axis     none
                scale is drawn if any.
                It should be one of
                "left", "right" or "none"

  -graph_height Specify height of the graph   Same as the
                                              "height" option.

=head1 BUGS

Please report them.

=head1 SEE ALSO

L<Bio::Graphics::Panel>,
L<Bio::Graphics::Track>,
L<Bio::Graphics::Glyph::transcript2>,
L<Bio::Graphics::Glyph::anchored_arrow>,
L<Bio::Graphics::Glyph::arrow>,
L<Bio::Graphics::Glyph::box>,
L<Bio::Graphics::Glyph::primers>,
L<Bio::Graphics::Glyph::segments>,
L<Bio::Graphics::Glyph::toomany>,
L<Bio::Graphics::Glyph::transcript>,

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

