package Bio::Graphics::Glyph::wiggle_box;

# $Id: wiggle_box.pm,v 1.2 2008/09/19 15:29:21 lstein Exp $

use strict;
use base qw(Bio::Graphics::Glyph::box Bio::Graphics::Glyph::smoothing);
use File::Spec;

sub draw {
  my $self = shift;
  my ($gd,$left,$top,$partno,$total_parts) = @_;
  my $feature   = $self->feature;

  my ($wigfile) = $feature->attributes('wigfile');
  if ($wigfile) {
    $self->draw_wigfile($feature,$self->rel2abs($wigfile),@_);
    $self->draw_label(@_)       if $self->option('label');
    $self->draw_description(@_) if $self->option('description');
    return;
  }

  my ($wigdata) = $feature->attributes('wigdata');
  if ($wigdata) {
    $self->draw_wigdata($feature,$wigdata,@_);
    $self->draw_label(@_)       if $self->option('label');
    $self->draw_description(@_) if $self->option('description');
    return;
  }

  return $self->SUPER::draw(@_);
}

sub wig {
  my $self = shift;
  my $d = $self->{wig};
  $self->{wig} = shift if @_;
  $d;
}

sub draw_wigdata {
    my $self    = shift;
    my $feature = shift;
    my $data    = shift;

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
    $self->_draw_wigfile($feature,$wig,@_);
}

sub draw_wigfile {
  my $self    = shift;
  my $feature = shift;
  my $wigfile = shift;

  eval "require Bio::Graphics::Wiggle" unless Bio::Graphics::Wiggle->can('new');
  my $wig = Bio::Graphics::Wiggle->new($wigfile) or die;
  $self->wig($wig);

  $self->_draw_wigfile($feature,$wig,@_);
}

sub _draw_wigfile {
    my $self    = shift;
    my $feature = shift;
    my $wig     = shift;
    my ($gd,$left,$top) = @_;

    my $start          = $self->smooth_start;
    my $end            = $self->smooth_end;

    my ($x1,$y1,$x2,$y2) = $self->bounds($left,$top);
    $self->draw_segment($gd,
			$start,$end,
			$wig,$start,$end,
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
  my $scale  = $self->scale;  # pixels per base pair
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

  return unless $data && ref $data && @$data > 0 && grep {$_} @$data;

  # allocate colors
  my $bg_idx = $self->panel->translate_color($self->panel->rgb($self->bgcolor));
  my $fg_idx = $self->panel->translate_color($self->panel->rgb($self->fgcolor)) || $bg_idx;

  $pixels_per_step = $scale * $step;
  $pixels_per_step = 1 if $pixels_per_step < 1;
  my $datapoints_per_base  = @$data/$length;
  my $pixels_per_datapoint = $self->panel->width/@$data * $data_width_ratio;

  my $xstart;
  for (my $i = 0; $i <= @$data ; $i++) {
    $xstart ||= $x1 + $pixels_per_datapoint * $i if $data->[$i];
    # trigger to draw the previous box is empty space of the end of the stack
    if (!$data->[$i] || ($i+1 == @$data)) {
      $xstart || next;
      my $xend = $x1 + $pixels_per_datapoint * $i;
      $self->filled_box($gd,$xstart,$y1,$xend,$y2,$bg_idx,$fg_idx);
      undef $xstart;
    }
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

Bio::Graphics::Glyph::wiggle_box - A generic box glyph compatible with dense "wig"data

=head1 SYNOPSIS

  See <Bio::Graphics::Panel> and <Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph works like the regular 'box' glyph but takes value data in
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

This glyph is intended for dense, qualitative feature data.  Any score data
for each data point is only evaluated for true/false, when true, a box
of the specified bgcolor is drawn, when false, nothing is drawn.  No
data smoothing is used.

Two primary benefits of using this glyph (with wiggle data) are:

 1) For large, genome-wide data sets, the speed of panel rendering is 
    greatly improved.
 2) Large sets of related features can be rendered as a UCSC-style subtrack
    without the need for  aggregation or a GFF3 containment hierarchy. 

A disadvantage to this approach is that individual features will have no
attributes associated with them and will appear as anonymous blocks within
a sub-track. 

An example use for this glyph is annotated transcribed regions from microarray
experiments.  Such regions are identified based on raw microarray data but do 
not necessarily have a score associated with them.  In this case, using the 
wiggle_box glyph provides a graphical summary of an expression array experiment.

=head2 DATA

The wiggle data used for this glyph should be loaded using the 'BED' format in
order to allow features of variable width.  The fourth column should be a true
value, with numeric or ".".  An example is shown below:

 track type=wiggle_0 name="transfrags" description="D. melanogaster transcribed fragments 0-2hrs"
 2L      9309    9451    1
 2L      10697   11021   1
 2L      11101   11345   1
 2L      11410   11521   1
 2L      11771   12243   1
 2L      12314   12954   1
 2L      13516   15746   1
 2L      17033   17191   1
 2L      18232   18580   1
 2L      19860   19999   1

=head2 OPTIONS

This glyph accepts the standard generic option set.  It differs in that
the label and description and title/mouseover labels apply to the whole, 
panel-wide sub-track feature rather than to individual boxes.

See Bio::Graphics::Glyph::wiggle_xyplot for a description of the
wiggle-specific options and data formats.

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

Sheldon McKay  E<lt>mckays@cshl.eduE<gt>.

Copyright (c) 2008 Cold Spring Harbor Laboratory

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut
