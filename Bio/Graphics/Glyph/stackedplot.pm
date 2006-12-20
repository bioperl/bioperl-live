package Bio::Graphics::Glyph::stackedplot;

use strict;
use base 'Bio::Graphics::Glyph::xyplot';
use Carp 'cluck';

use constant TOP_SPACING => 12;

sub width_needed {
  my $self = shift;
  my $column_width   = $self->column_width;
  my $column_spacing = $self->column_spacing;
  my $scale_width    = $self->scale_width;
  my $columns        = $self->data_series;
  return (@$columns-2) * $column_width + (@$columns-1)*$column_spacing + $scale_width;
}


sub pad_top {
  my $self    = shift;
  my $top  = $self->SUPER::pad_top;
  return $top + TOP_SPACING + $self->delegate_height;
}
sub pad_bottom {
  my $self = shift;
  my @labels  = $self->category_labels;
  return $self->SUPER::pad_bottom unless @labels;
  return $self->font('gdTinyFont')->height;
}

sub column_width    { shift->option('column_width')     || 4  }
sub column_spacing  { shift->option('column_spacing')   || 2  }
sub delegate_height { shift->option('delegate_height')  || 8  }
sub scale_width     { shift->option('scale_width')      || 20 }

sub pad_left {
  my $self = shift;
  my $pad          = $self->SUPER::pad_left;
  my $width_needed = ($self->width_needed - $self->width)/2;
  return $pad > $width_needed ? $pad : $width_needed;
}

sub pad_right {
   my $self = shift;
   my $pad          = $self->SUPER::pad_right;
   my $width_needed = ($self->width_needed - $self->width)/2;
   return $pad > $width_needed ? $pad : $width_needed;
}

sub maxdepth { 0 }

# this behaves more like the image glyph -- it draws a generic glyph, two diagonal lines, and then the
# plot underneath.
sub draw {
   my $self = shift;
   my $gd       = shift;
   my ($dx,$dy) = @_;
   my($x1,$y1,$x2,$y2) = $self->bounds($dx,$dy);

   my $top = $y1 - $self->pad_top;
   my $bottom = $y2;
   $self->filled_box($gd,$x1,$top,$x2,$top+6);

  my $width        = $self->width_needed;
  my $graph_top    = $y1;
  my $xmid         = ($x1+$x2) / 2;
  my $graph_left   = $xmid - $width/2;
  my $graph_right  = $xmid + $width/2;
  my $fgcolor      = $self->fgcolor;

  if (TOP_SPACING > 0) {
    $top += 6;
    $gd->line($x1,$top+2,$x1,$top+4,$fgcolor);
    $gd->line($x2,$top+2,$x2,$top+4,$fgcolor);
    $gd->line($x1,$top+4,$graph_left,$y1-4,$fgcolor);
    $gd->line($x2,$top+4,$graph_right,$y1-4,$fgcolor);
    $gd->line($graph_left,$y1-4,$graph_left,$y1-2,$fgcolor);
    $gd->line($graph_right,$y1-4,$graph_right,$y1-2,$fgcolor);
  }

   my $min_score = $self->option('min_score') || 0.0;
   my $max_score = $self->option('max_score') || 1.0;
   my $height = $self->height;
   my $scale  = $max_score > $min_score ? $height/($max_score-$min_score) : 1;
   my $y_origin = $min_score <= 0 ? $bottom - (0 - $min_score) * $scale : $bottom;
   $y_origin    = $top if $max_score < 0;
   $self->_draw_scale($gd,$scale,$min_score,$max_score,$dx+$self->pad_left,$dy,$y_origin);
   $self->draw_stackedplot($gd,$dx,$dy,$scale,$min_score,$max_score);
#   $self->SUPER::draw($gd,$dx,$dy);
}

sub draw_stackedplot {
  my $self = shift;
  my ($gd,$left,$top,$scale,$min,$max) = @_;

  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;
  my @colors  = $self->series_colors;
  my @labels  = $self->category_labels;
  my $column_width   = $self->column_width;
  my $column_spacing = $self->column_spacing;
  my $font      = $self->font('gdTinyFont');
  my $fwidth    = $font->width;
  my $fontcolor = $self->fontcolor;

  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries($left,$top);

  # data_series() returns 1 or more values to stack upwards
  # the totals of the values must be no greater than max_score
  if (my $values = $self->data_series) {
    my $x_offset = -$self->pad_left;
    $gd->line($x1+$x_offset,$y2,$x1+$self->pad_left,$y2,$fgcolor);

    for (my $cluster = 0; $cluster < @$values; $cluster++) {
      # this will give us a series of data series
      my $series = $values->[$cluster];
      my $y_offset = 0;
      for (my $i = 0; $i < @$series; $i++) {
	my $value = $series->[$i];
	my $v      = $self->clip($value,$min,$max);	
	my $color  = $colors[$i] || $bgcolor;

	my $y      = $y2 - ($v-$min) * $scale;
	my $box_bottom = $y2 - $y_offset;
	my $box_top    = $y  - $y_offset;
	$self->filled_box($gd,$x1+$x_offset,$box_top,$x1+$column_width+$x_offset,$box_bottom,$color);
	$y_offset += $box_bottom-$box_top;
      }
      if (@labels) {
	my $x = $x1+$x_offset+($column_width-$fwidth*$labels[$cluster])/2-1;
	$gd->string($font,$x,$y2,$labels[$cluster],$fontcolor);
      }
      $x_offset += $column_spacing+$column_width;
    }

  }
}

sub clip {
  my $self = shift;
  my ($value,$min,$max) = @_;
  $value = $min if defined $min && $value < $min;
  $value = $max if defined $max && $value > $max;
  return $value;
}

sub series_colors {
  my $self    = shift;
  my $values  = $self->option('series_colors');
  my @colors;
  if ($values && !ref $values) {
    @colors = split /\s+/,$values;
  } elsif (ref $values eq 'ARRAY') {
    @colors = @$values;
  } else {
    @colors = qw(red blue green orange brown grey black);
  }
  return map {$self->factory->translate_color($_)} @colors;
}

sub category_labels {
  my $self    = shift;
  my $values  = $self->option('category_labels');
  my @labels;
  if ($values && !ref $values) {
    @labels = split /\s+/,$values;
  } elsif (ref $values eq 'ARRAY') {
    @labels = @$values;
  }
  return @labels;
}

# NOTE!
# probably data_series should return this:
# [series1 => [value1,value2,value3,value4],
#  series2 => [value1,value2,value3,value4],
#  series3 => [value1,value2,value3,value4],
#  ...
# ]
sub data_series {
  my $self  = shift;
  my $values = $self->option('series');
  return $values if defined $values;

  # otherwise get it from the feature
  my @values;
  my @tagvalues = $self->feature->get_tag_values('series');
  for my $v (@tagvalues) {
    if (ref $v && ref $v eq 'ARRAY') {  # already in right format
      push @values,$v;
    } else {
      push @values,[split /[,\w]/,$v];
    }
  }
  return \@values;
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::stackedplot - The stackedplot glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

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

  -hilite       Highlight color                undef (no color)

In addition, the alignment glyph recognizes all the options of the
xyplot glyph, as well as the following glyph-specific option:

  Option         Description                  Default
  ------         -----------                  -------


=head1 EXAMPLES


=back

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

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

