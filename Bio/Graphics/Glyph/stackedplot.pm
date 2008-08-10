package Bio::Graphics::Glyph::stackedplot;

use strict;
use base 'Bio::Graphics::Glyph::fixedwidth';
use GD::Simple;
use Carp 'cluck';
use Memoize;
memoize('scale_width');
memoize('width_needed');

sub width_needed {
  my $self = shift;
  my $column_width   = $self->column_width;
  my $column_spacing = $self->column_spacing;
  my $columns        = $self->data_series;
  my $needed          = @$columns * ($column_width + $column_spacing);
  return $needed;
}

#sub pad_right {
#  my $self = shift;
#  my $pr   = $self->SUPER::pad_right;
#  my $sw   = $self->scale_width;
#  return $pr+$sw;
#}
sub pad_right {
  my $self = shift;
  my $pr      = $self->SUPER::pad_right;
  my $sw      = $self->scale_width + $self->column_width/2;
  my $content_width = $self->width_needed;
  my $total   = $sw + $content_width;
  my $additional = $pr - $total;
  return $pr if $pr > $additional;
  return $additional;
}

sub scale_width {
  my $self = shift;
  return 0 unless $self->do_draw_scale;
  my ($min,$max) = $self->min_max;
  my $middle     = ($min+$max)/2;
  my ($longest)  = sort {$b<=>$a} map {length($_)} ($min,$middle,$max);
  return $longest * $self->scale_font->width;
}

sub pad_bottom {
  my $self = shift;
  my @labels  = $self->column_labels;
  my $bottom = $self->SUPER::pad_bottom;
  return $bottom unless @labels;
  return $bottom + $self->font('gdTinyFont')->height;
}

sub column_width    { shift->option('column_width')     || 8  }
sub column_spacing  { shift->option('column_spacing')   || 2  }
sub maxdepth { 0 }

sub min_max {
  my $self = shift;
  my $min_score = $self->option('min_score') || 0.0;
  my $max_score = $self->option('max_score') || 1.0;
  return ($min_score,$max_score);
}

# this behaves more like the image glyph -- it draws a generic glyph, two diagonal lines, and then the
# plot underneath.
sub draw_contents {
   my $self = shift;
   my ($gd,$left,$top,$right,$bottom) = @_;

   my ($min_score,$max_score) = $self->min_max;
   my $height    = $bottom-$top;

   my $scale    = $max_score > $min_score ? $height/($max_score-$min_score) : 1;
   my $y_origin = $min_score <= 0 ? $bottom - (0 - $min_score) * $scale : $bottom;
   my $scale_width = $self->scale_width;
   $y_origin    = $top if $max_score < 0;
   $self->draw_stackedplot($gd,$left,$right,$top,$y_origin,$scale,$min_score,$max_score);

   $self->draw_scale($gd,$left,$top,$right,$bottom);
}

sub draw_stackedplot {
  my $self = shift;
  my ($gd,$left,$right,$top,$bottom,$scale,$min,$max) = @_;

  my $fgcolor = $self->fgcolor;
  my $bgcolor = $self->bgcolor;
  my @colors  = $self->series_colors;
  my @labels  = $self->column_labels;
  my $column_width   = $self->column_width;
  my $column_spacing = $self->column_spacing;
  my $font      = $self->column_font;
  my $fwidth    = $font->width;
  my $fontcolor = $self->fontcolor;

  # data_series() returns 1 or more values to stack upwards
  # the totals of the values must be no greater than max_score
  if (my $values = $self->data_series) {
    my $x_offset = 0;

    for (my $cluster = 0; $cluster < @$values; $cluster++) {
      # this will give us a series of data series
      my $series = $values->[$cluster];
      my $y_offset = 0;
      for (my $i = 0; $i < @$series; $i++) {
	my $value = $series->[$i];
	my $v      = $self->clip($value,$min,$max);
	my $color  = $colors[$i] || $bgcolor;

	my $y          = $bottom - ($v-$min) * $scale;
	my $box_bottom = $bottom - $y_offset;
	my $box_top    = $y      - $y_offset;
	$self->filled_box($gd,$left+$x_offset,$box_top,$left+$column_width+$x_offset,$box_bottom,$color);
	$y_offset += $box_bottom-$box_top;
      }
      if (@labels) {
	my $x = $left+$x_offset+($column_width-$fwidth*$labels[$cluster])/2-1;
	$gd->string($font,$x,$bottom,$labels[$cluster],$fontcolor);
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
  return map {$self->translate_color($_)} @colors;
}

sub column_labels {
  my $self    = shift;
  my $values  = $self->option('column_labels');
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
      push @values,[split /[,\s]+/,$v];
    }
  }
  return \@values;
}
sub do_draw_scale {
  my $self = shift;
  my $drawit = $self->option('draw_scale');
  return defined $drawit ? $drawit : 1;
}
sub scale_font {
  my $self = shift;
  $self->getfont('scale_font','gdTinyFont');
}
sub column_font {
  my $self = shift;
  $self->getfont('column_font','gdSmallFont');
}

sub draw_scale {
  my $self = shift;
  my ($gd,$left,$top,$right,$bottom) = @_;
  my ($min,$max) = $self->min_max;

  $self->panel->startGroup($gd);

  my $simple = GD::Simple->new($gd);
  $simple->font($self->scale_font);
  my $dx     = 1;
  my $dy     = $simple->font->height/2;

  # these drew a vertical scale line, which didn't look very nice
#  $simple->moveTo($right,$bottom);
#  $simple->lineTo($right,$top);

  $simple->moveTo($right,$top);
  $simple->line(3);
  $simple->move($dx,$dy);
  $simple->string($max);

  $simple->moveTo($right,($bottom+$top)/2);
  $simple->line(3);
  $simple->move($dx,$dy);
  $simple->string(($max+$min)/2);

  $simple->moveTo($right,$bottom);
  $simple->line(3);
  $simple->move($dx,$dy);
  $simple->string($min);

  $self->panel->endGroup($gd);
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::stackedplot - The stackedplot glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

The stackedplot glyph can be used to draw quantitative feature data
using a stacked column plot. It differs from the xyplot glyph in that
the plot applies to a single top level feature, not a group of
subfeatures. The data to be graphed is derived from an attribute
called "data_series."

The data to be graphed is represented as a list of arrays:

 (
 [1, 2, 8],
 [6, 1, 1],
 [10,8, 0],
 [1, 1, 1],
 )

Each array is a column in the stacked plot. Its values become the
subdivisions of the column. In this example, there are four columns,
each of which has three subdivisions.

You can add labels to the columns and change the colors of the
subdivisions.

To assign data to a feature, you can add a "series" tag:

 $snp1    = Bio::SeqFeature::Generic ->new (-start     => 500,-end=>501,
			                    -display_name =>'example',
					    -tag=> { series => [
 							     [10,20,30],
 							     [30,30,0],
 							     [5,45,10],
 							     [5,45,10],
 							     [5,45,10],
 							     [50,0,50],
 							    ],
						     }
					       );

Note that the series tag must consist of an array of arrays.

If you are using a gff3 representation, you can load a database with
data that looks like this:

 chr1 test feature 1 1000 . . . series=10 20 30;series=30 30 0;series=5 45 10...

If you are using a gff2 representation, you can load a database with
data that looks like this:

 chr1 test feature 1 1000 . . . series 10 20 30; series 30 30 0 series 5 45 10...

Or you can pass a callback to the -series option:

 $panel->add_track(\@data,
		  -glyph     => 'stackedplot',
                  -series       => sub {
                                  my $feature = shift;
                             return [
			        [10,20,30],
 				[30,30,0],
 				[5,45,10],
                            ]
                          }
		 );

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

 -fixed_gap     Vertical distance between      8
                the rectangle that shows
                the start:end range of
                the feature and the fixed
                width stacked plot.

 -series_colors A list giving a series of     red,blue,green,orange,
                color names for the data      brown,grey,black
                series (the values inside
                each stacked column).

 -column_labels A list of labels to print     -none-
                underneath each column.

 -column_width  The width of each column.     8

 -column_spacing Spacing between each         2
                column.

 -min_score     Minimum score for the         0.0
                sum of the members of
                each data series.

 -max_score     Maximum score for the         1.0
                sum of the members of each
                data series.

 -scale_font    Font to use for the scale.    gdTinyFont

 -column_font   Font to use for the column    gdSmallFont
                labels.

 -draw_scale    Whether to draw a scale to    true
                right of the columns.

Note that -min_score and -max_score represent the minimum and maximum
SUM of all the values in the data series. For example, if your largest
column contains the series (10,20,30), then the -max_score is 60.

=head1 EXAMPLE

To understand how this glyph works, try running and modifying the following example:

 #!/usr/bin/perl

 use strict;
 use warnings;

 use Bio::Graphics;
 use Bio::SeqFeature::Generic;

 my $segment  = Bio::Graphics::Feature->new(-start=>1,-end=>700);

 my $snp1     = Bio::SeqFeature::Generic ->new (-start     => 500,-end=>590,
 					        -display_name =>'fred',
					        -tag=> { series => [
 								     [10,20,30],
 								     [30,30,0],
 								     [5,45,10],
 								     [5,45,10],
 								     [5,45,10],
 								     [50,0,50],
 								    ],
						      },
					        -source=>'A test',
					        );

 my $snp2     = Bio::SeqFeature::Generic->new(-start     => 300,
					      -end       => 301,
					      -display_name  => 'rs12345',
					      -tag=> {
						     series => [
								     [30,20,10 ],
								     [80,10,10 ],
							       ],
						    },
					      -source=>'Another test',
					    );

 my $panel = Bio::Graphics::Panel->new(-segment=>$segment,-width=>800);

 $panel->add_track($segment,-glyph=>'arrow',-double=>1,-tick=>2);
 $panel->add_track([$snp1,$snp2],
		   -height    => 50,
		   -glyph     => 'stackedplot',
		   -fixed_gap => 12,
		   -series_colors    => [qw(red blue lavender)],
		   -column_labels => [qw(a b c d e f g)],
		   -min_score => 0,
		   -max_score => 100,
		   -column_width => 8,
		   -column_font  => 'gdMediumBoldFont',
		   -scale_font   => 'gdTinyFont',
		   -label    => 1,
		   -description=>1,
		  );
 print $panel->png;

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

