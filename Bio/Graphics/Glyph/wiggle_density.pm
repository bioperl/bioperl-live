package Bio::Graphics::Glyph::wiggle_density;
# $Id: wiggle_density.pm,v 1.4 2008/09/24 23:34:43 lstein Exp $

use strict;
use base qw(Bio::Graphics::Glyph::box Bio::Graphics::Glyph::smoothing Bio::Graphics::Glyph::minmax);
use File::Spec;

sub draw {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $feature   = $self->feature;

  $self->panel->startGroup($gd);

  my ($wigfile) = $feature->attributes('wigfile');
  if ($wigfile) {
    $self->draw_wigfile($self->rel2abs($wigfile),@_);
    $self->draw_label(@_)       if $self->option('label');
    $self->draw_description(@_) if $self->option('description');
    return;
  }

  my ($wigdata) = $feature->attributes('wigdata');
  if ($wigdata) {
    $self->draw_wigdata($wigdata,@_);
    $self->draw_label(@_)       if $self->option('label');
    $self->draw_description(@_) if $self->option('description');
    return;
  }

  my ($densefile) = $feature->attributes('densefile');
  if ($densefile) {
    $self->draw_densefile($self->rel2abs($feature),$densefile,@_);
    $self->draw_label(@_)       if $self->option('label');
    $self->draw_description(@_) if $self->option('description');
    return;
  }

  $self->panel->endGroup($gd);

  return $self->SUPER::draw(@_);
}

sub wig {
  my $self = shift;
  my $d = $self->{wig};
  $self->{wig} = shift if @_;
  $d;
}

sub draw_wigfile {
  my $self    = shift;
  my $wigfile = shift;

  warn "rendering $wigfile";

  eval "require Bio::Graphics::Wiggle" unless Bio::Graphics::Wiggle->can('new');
  my $wig = eval { Bio::Graphics::Wiggle->new($wigfile)};

  unless ($wig) {
      warn $@;
      return $self->SUPER::draw(@_);
  }
  $self->wig($wig);

  $self->_draw_wigfile(@_);
}

sub draw_wigdata {
    my $self = shift;
    my $data = shift;

    eval "require MIME::Base64" 
	unless MIME::Base64->can('decode_base64');
    my $unencoded_data = MIME::Base64::decode_base64($data);

    my $wig = eval { Bio::Graphics::Wiggle->new() };
    unless ($wig) {
	warn $@;
	return $self->SUPER::draw(@_);
    }

    $wig->import_from_wif($unencoded_data);

    $self->wig($wig);
    $self->_draw_wigfile(@_);
}

sub _draw_wigfile {
    my $self = shift;
    my $wig  = $self->wig;
    my ($gd,$left,$top) = @_;

    my $smoothing      = $self->get_smoothing;
    my $smooth_window  = $self->smooth_window;
    my $start          = $self->smooth_start;
    my $end            = $self->smooth_end;

    $wig->window($smooth_window);
    $wig->smoothing($smoothing);
    my ($x1,$y1,$x2,$y2) = $self->bounds($left,$top);
    $self->draw_segment($gd,
			$start,$end,
			$wig,$start,$end,
			1,1,
			$x1,$y1,$x2,$y2);
}

sub draw_densefile {
  my $self = shift;
  my $feature   = shift;
  my $densefile = shift;
  my ($gd,$left,$top) = @_;

  my ($denseoffset) = $feature->attributes('denseoffset');
  my ($densesize)   = $feature->attributes('densesize');
  $denseoffset ||= 0;
  $densesize   ||= 1;

  my $smoothing      = $self->get_smoothing;
  my $smooth_window  = $self->smooth_window;
  my $start          = $self->smooth_start;
  my $end            = $self->smooth_end;

  my $fh         = IO::File->new($densefile) or die "can't open $densefile: $!";
  eval "require Bio::Graphics::DenseFeature" unless Bio::Graphics::DenseFeature->can('new');

  my $dense = Bio::Graphics::DenseFeature->new(-fh=>$fh,
					       -fh_offset => $denseoffset,
					       -start     => $feature->start,
					       -smooth    => $smoothing,
					       -recsize   => $densesize,
					       -window    => $smooth_window,
					      ) or die "Can't initialize DenseFeature: $!";

  my ($x1,$y1,$x2,$y2) = $self->bounds($left,$top);
  $self->draw_segment($gd,
		      $start,$end,
		      $dense,$start,$end,
		      1,1,
		      $x1,$y1,$x2,$y2);
}

sub draw_segment {
  my $self = shift;
  my ($gd,
      $start,$end,
      $seg_data,
      $seg_start,$seg_end,
      $step,$span,
      $x1,$y1,$x2,$y2) = @_;

  # clip, because wig files do no clipping
  $seg_start = $start      if $seg_start < $start;
  $seg_end   = $end        if $seg_end   > $end;

  # figure out where we're going to start
  my $scale           = $self->scale;  # pixels per base pair
  my $pixels_per_span = $scale * $span + 1;
  my $pixels_per_step = 1;
  my $length          = $end-$start+1;

  # if the feature starts before the data starts, then we need to draw
  # a line indicating missing data (this only happens if something went
  # wrong upstream)
  if ($seg_start > $start) {
    my $terminus = $self->map_pt($seg_start);
    $start = $seg_start;
    $x1    = $terminus;
  }
  # if the data ends before the feature ends, then we need to draw
  # a line indicating missing data (this only happens if something went
  # wrong upstream)
  if ($seg_end < $end) {
    my $terminus = $self->map_pt($seg_end);
    $end = $seg_end;
    $x2    = $terminus;
  }

  return unless $start < $end;

  # get data values across the area
  my $samples = $length < $self->panel->width ? $length : $self->panel->width;
  my $data    = $seg_data->values($start,$end,$samples);

  # scale the glyph if the data end before the panel does
  my $data_width = $end - $start;
  my $data_width_ratio;
  if ($data_width < $self->panel->length) {
    $data_width_ratio = $data_width/$self->panel->length;
  }
  else {
    $data_width_ratio = 1;
  }

  return unless $data && ref $data && @$data > 0;


  my $min_value = $self->min_score;
  my $max_value = $self->max_score;

  unless (defined $min_value && defined $max_value) {
      my ($min,$max) = $self->minmax($data);
      $min_value ||= $min;
      $max_value ||= $max;
  }

  # allocate colors
  # There are two ways to do this. One is a scale from min to max. The other is a
  # bipartite scale using one color range from zero to min, and another color range
  # from 0 to max. The latter behavior is triggered when the config file contains
  # entries for "pos_color" and "neg_color" and the data ranges from < 0 to > 0.

  my $poscolor       = $self->pos_color;
  my $negcolor       = $self->neg_color;

  my $data_midpoint  =   $self->bicolor_pivot;
  my $bicolor   = $poscolor != $negcolor
                       && $min_value < $data_midpoint
		       && $max_value > $data_midpoint;

  my ($rgb_pos,$rgb_neg,$rgb);
  if ($bicolor) {
      $rgb_pos = [$self->panel->rgb($poscolor)];
      $rgb_neg = [$self->panel->rgb($negcolor)];
  } else {
      $rgb = [$self->panel->rgb($self->bgcolor)];
  }


  my %color_cache;

  @$data = reverse @$data if $self->flip;

  if (@$data <= $self->panel->width) { # data fits in width, so just draw it

    $pixels_per_step = $scale * $step;
    $pixels_per_step = 1 if $pixels_per_step < 1;
    my $datapoints_per_base  = @$data/$length;
    my $pixels_per_datapoint = $self->panel->width/@$data * $data_width_ratio;

    for (my $i = 0; $i <= @$data ; $i++) {
      my $x          = $x1 + $pixels_per_datapoint * $i;
      my $data_point = $data->[$i] || next;
      $data_point    = $min_value if $min_value > $data_point;
      $data_point    = $max_value if $max_value < $data_point;
      my ($r,$g,$b)  = $bicolor
	  ? $data_point > $data_midpoint ? $self->calculate_color($data_point,$rgb_pos,
								  $data_midpoint,$max_value)
	                                 : $self->calculate_color($data_point,$rgb_neg,
								  $data_midpoint,$min_value)
          : $self->calculate_color($data_point,$rgb,
				   $min_value,$max_value);
      my $idx        = $color_cache{$r,$g,$b} ||= $self->panel->translate_color($r,$g,$b);
      $self->filled_box($gd,$x,$y1,$x+$pixels_per_datapoint,$y2,$idx,$idx);
    }

  } else {     # use Sheldon's code to subsample data
      $pixels_per_step = $scale * $step;
      my $pixels = 0;

      # only draw boxes 2 pixels wide, so take the mean value
      # for n data points that span a 2 pixel interval
      my $binsize = 2/$pixels_per_step;
      my $pixelstep = $pixels_per_step;
      $pixels_per_step *= $binsize;
      $pixels_per_step *= $data_width_ratio;
      $pixels_per_span = 2;

      my $scores = 0;
      my $defined;

      for (my $i = $start; $i < $end ; $i += $step) {
	# draw the box if we have accumulated >= 2 pixel's worth of data.
	if ($pixels >= 2) {
	  my $data_point = $defined ? $scores/$defined : 0;
	  $scores  = 0;
	  $defined = 0;

	  $data_point    = $min_value if $min_value > $data_point;
	  $data_point    = $max_value if $max_value < $data_point;
	  my ($r,$g,$b)  = $bicolor
	      ? $data_point > $data_midpoint ? $self->calculate_color($data_point,$rgb_pos,
								      $data_midpoint,$max_value)
	                                     : $self->calculate_color($data_point,$rgb_neg,
								      $data_midpoint,$min_value)
	      : $self->calculate_color($data_point,$rgb,
				       $min_value,$max_value);
	  my $idx        = $color_cache{$r,$g,$b} ||= $self->panel->translate_color($r,$g,$b);
	  $self->filled_box($gd,$x1,$y1,$x1+$pixels_per_span,$y2,$idx,$idx);
	  $x1 += $pixels;
	  $pixels = 0;
	}

	my $val = shift @$data;
	# don't include undef scores in the mean calculation
	# $scores is the numerator; $defined is the denominator
	$scores += $val if defined $val;
	$defined++ if defined $val;

	# keep incrementing until we exceed 2 pixels
	# the step is a fraction of a pixel, not an integer
	$pixels += $pixelstep;
      }
  }
}

sub series_mean {
    my $self = shift;
    my $wig = $self->wig or return;
    return eval {$wig->mean} || undef;
}

sub calculate_color {
  my $self = shift;
  my ($s,$rgb,$min_score,$max_score) = @_;
  $s ||= $min_score;

  my $relative_score = ($s-$min_score)/($max_score-$min_score);
  $relative_score -= .1 if $relative_score == 1;
  return map { int(254.9 - (255-$_) * min(max( $relative_score, 0), 1)) } @$rgb;
}

sub min { $_[0] < $_[1] ? $_[0] : $_[1] }
sub max { $_[0] > $_[1] ? $_[0] : $_[1] }

sub minmax {
  my $self = shift;
  my $data = shift;

  if (my $wig = $self->wig) {
    return ($wig->min,$wig->max);
  } else {
      return $self->SUPER::minmax($data);
  }
}

sub rel2abs {
    my $self = shift;
    my $wig  = shift;
    my $path = $self->option('basedir');
    return File::Spec->rel2abs($wig,$path);
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::wiggle_density - A density plot compatible with dense "wig"data

=head1 SYNOPSIS

  See <Bio::Graphics::Panel> and <Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph works like the regular density but takes value data in
Bio::Graphics::Wiggle file format:

 reference = chr1
 ChipCHIP Feature1 1..10000 wigfile=./test.wig;wigstart=0
 ChipCHIP Feature2 10001..20000 wigfile=./test.wig;wigstart=656
 ChipCHIP Feature3 25001..35000 wigfile=./test.wig;wigstart=1312

The "wigfile" attribute gives a relative or absolute pathname to a
Bio::Graphics::Wiggle format file. The optional "wigstart" option
gives the offset to the start of the data. If not specified, a linear
search will be used to find the data. The data consist of a packed
binary representation of the values in the feature, using a constant
step such as present in tiling array data.

=head2 OPTIONS

The same as the regular graded_segments glyph, except that the
following options are recognized:

   Name        Value        Description
   ----        -----        -----------

   basedir     path         Path to be used to resolve "wigfile" and "densefile"
                                tags giving relative paths. Default is to use the
                                current working directory. Absolute wigfile &
                                densefile paths will not be changed.

   smoothing   method name  Smoothing method: one of "mean", "max", "min" or "none"

   smoothing_window 
               integer      Number of values across which data should be smoothed.

   bicolor_pivot
               name         Where to pivot the two colors when drawing bicolor plots.
                               Options are "mean" and "zero". A numeric value can
                               also be provided.

   pos_color   color        When drawing bicolor plots, the fill color to use for values
                              that are above the pivot point.

   neg_color   color        When drawing bicolor plots, the fill color to use for values
                              that are below the pivot point.

=head2 SPECIAL FEATURE TAGS

The glyph expects one or more of the following tags (attributes) in
feature it renders:

   Name        Value        Description
   ----        -----        -----------

   wigfile     path name    Path to the Bio::Graphics::Wiggle file for vales.
                            (required)

   densefile   path name    Path to a Bio::Graphics::DenseFeature object
                               (deprecated)

   denseoffset integer      Integer offset to where the data begins in the
                               Bio::Graphics::DenseFeature file (deprecated)

   densesize   integer      Integer size of the data in the Bio::Graphics::DenseFeature
                               file (deprecated)

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
L<Bio::Graphics::Glyph::allele_tower>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>steinl@cshl.eduE<gt>.

Copyright (c) 2007 Cold Spring Harbor Laboratory

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut
