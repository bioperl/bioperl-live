package Bio::Graphics::Glyph::extending_arrow;
# package to use for drawing an arrow, but use dash line
# to indicate feature goes beyond canvas

# *** need to figure out how to return glyph actual drawnout coords for panel->boxes

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

sub box {
  my $self = shift;
  my ($x1, $y1, $x2, $y2) = ($self->left, $self->top, $self->right, $self->bottom);
  $x1 = $self->panel->left  if $x1 < $self->panel->left;

  #figure out how much it was extended (except dx)
  my $start    = $self->feature->start;
  my $stop     = $start + $self->feature->length;
  my $p_seg_start = $self->panel->offset;
  my $pl =  $self->panel->pad_left;
  if ($start < $p_seg_start && $x1 eq $self->panel->left) {
      $x1 -= $self->panel->pad_left;
      my $map_x1 = $pl - $self->scale * ($p_seg_start - $start);
      if ($map_x1 > $x1){
          $x1 = $map_x1;
      }
  }
  if ($stop > $p_seg_start + $self->panel->length) {
      $x2 += $self->panel->pad_right;
      my $map_x2 = $pl + $self->scale * ($stop - $p_seg_start);
      if ($map_x2 <= $x2) {
          $x2 = $map_x2;
      }
  }
  return ($x1, $y1, $x2, $y2);
}

sub extended_left {
  my $self = shift;
  my $x1 = shift;
  my $dx = shift || 0;

  $x1 = $self->panel->left  if $x1 < $self->panel->left;
  #figure out if to extending and how much
  my $start    = $self->feature->start;
  my $stop     = $start + $self->feature->length;
  my $p_seg_start = $self->panel->offset;

  my $pl =  $self->panel->pad_left;
  if ($start < $p_seg_start && $x1 eq $self->panel->left) {
      $x1 = $pl + $dx - $self->scale * ($p_seg_start - $start);
  }
  return $x1;
}
sub extended_right {
  my $self = shift;
  my $x2 = shift;
  my $dx = shift || 0;
  $x2 = $self->panel->right if $x2 > $self->panel->right;
  #figure out if to extending and how much
  my $start    = $self->feature->start;
  my $stop     = $start + $self->feature->length;
  my $p_seg_start = $self->panel->offset;

  my $pl =  $self->panel->pad_left;
  if ($stop > $p_seg_start + $self->panel->length) {
      $x2 = $pl + $dx + $self->scale * ($stop - $p_seg_start);
  }
  return $x2;
}

sub dashline_length {
  my $self = shift;
  return $self->option('dashline_length') || 20;
}

sub draw_label {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $label = $self->label or return;

  my ($x1,$y1,$x2,$y2) = $self->bounds($left, $top);

  $x1 = $self->panel->left  if $x1 < $self->panel->left;

  my $extended_x1 = $self->extended_left($x1, $left);
  if ($extended_x1 < $x1) {
      $x1 -= $self->panel->pad_left;
      if ($extended_x1 > $x1) {
          $x1 = $extended_x1;
      }
  }
#  if ($start < $p_seg_start && $x1 eq $self->panel->left) {
#      my $pl =  $self->panel->pad_left;
#      $x1 -= $self->panel->pad_left;
#      my $map_x1 = $pl + $left - $self->scale * ($p_seg_start - $start);
#      $x1 = $map_x1 if ($map_x1 > $x1);
#  }
  my $font = $self->option('labelfont') || $self->font;
  $gd->string($font,
              $x1,
              $self->top + $top,
              $label,
              $self->fontcolor);
}

# override draw method
sub draw {
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

  my $dash_len = $self->dashline_length;
  my ($line_x1, $line_x2) = ($x1, $x2);
  my $extended_x1 = $self->extended_left($x1, $dx);
  if ($extended_x1 < $x1) {
      $x1 -= $self->panel->pad_left;
      if ($extended_x1 > $x1) {
          $x1 = $extended_x1;
          $line_x1 = $x1;
      } else {
          my $dash_x2 = $x1 + $dash_len;
          $dash_x2 += $a2 if $sw;
          $dash_x2 = $line_x2 if ($dash_x2 > $x2);
          $gd->dashedLine($x1,$center,$dash_x2,$center,$fg);
          $line_x1 = $dash_x2;
      }
  }
  my $extending_x2 = $self->extended_right($x2, $dx);
  if ($extending_x2 > $x2) {
      $x2 += $self->panel->pad_right - $a2;
      if ($extending_x2 <= $x2) {
          $x2 = $extending_x2;
          $line_x2 = $x2;
      } else {
          my $dash_x1 = $x2 - $dash_len;
          $dash_x1 -= $a2 if $ne;
          $dash_x1 = $x1 if ($dash_x1 < $x1);
          $gd->dashedLine($dash_x1,$center,$x2,$center,$fg);
          $line_x2 = $dash_x1;
      }
  }

  $gd->line($line_x1,$center,$line_x2,$center,$fg);
  $self->arrowhead($gd,$x1,$center,$a2,-1) if $sw; # west arrow
  $self->arrowhead($gd,$x2,$center,$a2,+1) if $ne; # east arrow
  $gd->line($x2,$center-$a2,$x2,$center+$a2,$fg) if $base_e; #east base
  $gd->line($x1,$center-$a2,$x1,$center+$a2,$fg) if $base_w; #west base

  # turn on ticks
  if ($self->option('tick')) {
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
      $gd->string($font,$middle,$center+$a2-1,$label,$font_color);
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
  }
  return ($sw,$ne,0,0) unless $self->option('base');
  return ($sw,$ne,!$sw,!$ne);
}

1;

=head1 NAME

Bio::Graphics::Glyph::extending_arrow -- The "extending arrow" glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph was designed to show a segment that goes beyond the panel.
Dashed line indicates the end goes beyond the panel and arrow
indicates the direction.

Also see the anchored_arrw and arrow glyphs.

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

In addition to the generic options, this glyph recognizes:

 Option Name         Description                 Default
 -----------         -----------                 -------

 dashline_length     length of drawn dash line   20

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

Shengqiang Shu

Copyright (c) 2001 Berkeley Drosophila Genome Project

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
