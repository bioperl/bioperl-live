package Bio::Graphics::Glyph::ruler_arrow;
# package to use for drawing an arrow as ruler (5' and 3' are marked as label)

use strict;
use vars '@ISA';
use Bio::Graphics::Glyph::generic;
@ISA = 'Bio::Graphics::Glyph::generic';

my %UNITS = (K => 1000,
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
  $self->draw_label(@_) if ($self->option('label'));
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
#  $self->draw_label($gd,$dx,$dy) if ($self->option('label') && !$self->option('ruler'));
  # this draws the label aligned to the left
}

sub draw_parallel {
  my $self = shift;
  my $gd = shift;
  my ($dx,$dy) = @_;
  my ($x1,$y1,$x2,$y2) = $self->bounds(@_);

  my $fg = $self->set_pen;
  my $a2 = ($self->height)/2;
  my $center = $y1+$a2;

  $x1 = $self->panel->left  if $x1 < $self->panel->left;
  $x2 = $self->panel->right if $x2 > $self->panel->right;

  my ($sw,$ne,$base_w,$base_e) = $self->arrowheads;
  $gd->line($x1,$center,$x2,$center,$fg);
  $self->arrowhead($gd,$x1,$center,$a2,-1) if $sw; # west arrow
  $self->arrowhead($gd,$x2,$center,$a2,+1) if $ne; # east arrow
  $gd->line($x2,$center-$a2,$x2,$center+$a2,$fg) if $base_e; #east base
  $gd->line($x1,$center-$a2,$x1,$center+$a2,$fg) if $base_w; #west base

  # turn on ticks
  if ($self->option('tick')) {
      local $^W = 0;  # dumb uninitialized variable warning
    my $font = $self->font;
    my $width      = $font->width;
    my $font_color = $self->fontcolor;
    my $height   = $self->height;

    my $relative = $self->option('relative_coords');
    my $start    = $relative ? 1 : $self->feature->start;
    my $stop     = $start + $self->feature->length  - 1;

    my $offset   = $relative ? $self->feature->start-1 : 0;
    my $reversed = $self->feature->strand < 0;

    my $units = $self->option('units') || '';
    my $divisor = $UNITS{$units} || 1 if $units;

    my ($major_ticks,$minor_ticks) = $self->panel->ticks($start,$stop,$font,$divisor);

    ## Does the user want to override the internal scale?
    my $scale = $self->option('scale');

    my $left  = $sw ? $x1+$height : $x1;
    my $right = $ne ? $x2-$height : $x2;

    my $format = ($major_ticks->[1]-$major_ticks->[0])/($divisor||1) < 1 ? "%.1f$units" : "%d$units";

    for my $i (@$major_ticks) {
      my $tickpos = $dx + ($reversed ? $self->map_pt($stop - $i + $offset)
	                             : $self->map_pt($i + $offset));
      next if $tickpos < $left or $tickpos > $right;
      $gd->line($tickpos,$center-$a2,$tickpos,$center+$a2,$fg);
      my $label = $scale ? $i / $scale : $i;
      if ($units) {
	my $scaled = $label/$divisor;
	$label = sprintf($format,$scaled);
      }
      my $middle = $tickpos - (length($label) * $width)/2;
      $gd->string($font,$middle,$center+$a2-1,$label,$font_color)
        unless ($self->option('no_tick_label'));
    }

    if ($self->option('tick') >= 2) {
      my $a4 = $self->height/4;
      for my $i (@$minor_ticks) {
	my $tickpos = $dx + ($reversed ? $self->map_pt($stop - $i + $offset)
	                               : $self->map_pt($i + $offset));
	next if $tickpos < $left or $tickpos > $right;
	$gd->line($tickpos,$center-$a4,$tickpos,$center+$a4,$fg);
      }
    }
  }

  # add a label if requested
#  $self->draw_label($gd,$dx,$dy)       if ($self->option('label');
#  $self->draw_description($gd,$dx,$dy) if $self->option('description');
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
  }
  return ($sw,$ne,0,0) unless $self->option('base');
  return ($sw,$ne,!$sw,!$ne);
}

sub draw_label {
  my $self = shift;
  my ($gd,$left,$top) = @_;

  my $label5 = "5'";
  my $label3 = "3'";
  my $relative = $self->option('relative_coords');
  my $start    = $relative ? 1 : $self->feature->start;
  my $stop     = $start + $self->feature->length  - 1;

  my $offset   = $relative ? $self->feature->start-1 : 0;
  my $reversed = $self->feature->strand < 0;

  my $units = $self->option('units') || '';
  my $divisor = $UNITS{$units} || 1 if $units;

  my ($major_ticks,$minor_ticks) = $self->panel->ticks($start,$stop,$self->font,$divisor);
  my $tick_scale = " ($major_ticks bp/";
  $tick_scale .= ($self->option('tick') >= 2)?"major tick)":"tick)";

  my $top_left_label = $label5;
  $top_left_label .= $tick_scale if ($self->option('no_tick_label') && $self->option('tick'));
  #-1 direction mean lower end is 3' (minus strand on top)
  ($label5, $label3) = ($label3, $label5) if ($self->option('direction') == -1);
  my $x = $self->left + $left;
  $x = $self->panel->left + 1 if $x <= $self->panel->left;
  my $font = $self->option('labelfont') || $self->font;
  $gd->string($font,
              $x,
              $self->top + $top,
              $top_left_label,
              $self->fontcolor);
  my $x1 = $left + $self->panel->right - $font->width*length($label3);
  $gd->string($font,
              $x1,
              $self->top + $top,
              $label3,
              $self->fontcolor);
  if ($self->option('both')) {#minus strand as well
      $gd->string($font,
                  $x,
                  $self->bottom - $self->pad_bottom + $top,
                  $label3,
                  $self->fontcolor);
      my $x1 = $left + $self->panel->right - $font->width*length($label5);
      $gd->string($font,
                  $x1,
                  $self->bottom - $self->pad_bottom + $top,
                  $label5,
                  $self->fontcolor);
  }
}


1;


__END__

=head1 NAME

Bio::Graphics::Glyph::arrow - The "ruler_arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph draws arrows.  Label, if requested, will be 5' and 3' at both ends
and tick scale is printed if no_tick_label option is set and tick option set.
Depending on options, the arrows can be labeled, be oriented vertically 
or horizontally, or can contain major and minor ticks suitable for use as a scale.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description               Default
  ------      -----------               -------

  -tick       Whether to draw major         0
              and minor ticks.
	      0 = no ticks
	      1 = major ticks
	      2 = minor ticks
  -label      5' at start, 3' at end        0
              above arrow
  -both       5', 3' above,                 0
              and 3', 5' below arrow
  -direction  0 = ruler is plus strand      0
              -1 = ruler is minus strand

  -parallel   Whether to draw the arrow     true
	      parallel to the sequence
	      or perpendicular to it.

  -northeast  Force a north or east         true
	      arrowhead(depending 
	      on orientation)

  -east       synonym of above

  -southwest  Force a south or west         true
	      arrowhead(depending 
	      on orientation)

  -west       synonym of above

  -double     force-doubleheaded arrow

  -base       Draw a vertical base at the   false
              non-arrowhead side

  -scale      Reset the labels on the arrow false
              to reflect an externally 
              established scale.

Set -parallel to false to display a point-like feature such as a
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

Shengqiang Shu E<lt>sshu@bdgp.lbl.govE<gt>
Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 BDGP, Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
