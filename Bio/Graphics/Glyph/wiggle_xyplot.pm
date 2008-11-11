package Bio::Graphics::Glyph::wiggle_xyplot;

use strict;
use base qw(Bio::Graphics::Glyph::xyplot Bio::Graphics::Glyph::smoothing);
use IO::File;
use File::Spec;

# we override the draw method so that it dynamically creates the parts needed
# from the wig file rather than trying to fetch them from the database
sub draw {
  my $self = shift;
  my ($gd,$dx,$dy) = @_;

  my $feature     = $self->feature;
  my ($wigfile)   = $feature->attributes('wigfile');
  return $self->draw_wigfile($feature,$self->rel2abs($wigfile),@_) if $wigfile;

  my ($wigdata) = $feature->attributes('wigdata');
  return $self->draw_wigdata($feature,$wigdata,@_) if $wigdata;

  my ($densefile) = $feature->attributes('densefile');
  return $self->draw_densefile($feature,$self->rel2abs($densefile),@_) if $densefile;

  return $self->SUPER::draw(@_);
}

sub draw_wigfile {
  my $self = shift;
  my $feature = shift;
  my $wigfile = shift;

  eval "require Bio::Graphics::Wiggle" unless Bio::Graphics::Wiggle->can('new');
  my $wig = eval { Bio::Graphics::Wiggle->new($wigfile) };
  unless ($wig) {
      warn $@;
      return $self->SUPER::draw(@_);
  }
  $self->_draw_wigfile($feature,$wig,@_);
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

    $self->_draw_wigfile($feature,$wig,@_);
}

sub _draw_wigfile {
    my $self    = shift;
    my $feature = shift;
    my $wig     = shift;

    $wig->smoothing($self->get_smoothing);
    $wig->window($self->smooth_window);

    my $panel_start = $self->panel->start;
    my $panel_end   = $self->panel->end;
    my $start       = $feature->start > $panel_start ? $feature->start : $panel_start;
    my $end         = $feature->end   < $panel_end   ? $feature->end   : $panel_end;

    $self->wig($wig);
    $self->create_parts_for_dense_feature($wig,$start,$end);

    $self->SUPER::draw(@_);
}

sub wig {
  my $self = shift;
  my $d = $self->{wig};
  $self->{wig} = shift if @_;
  $d;
}

sub series_mean {
    my $self = shift;
    my $wig = $self->wig or return;
    return eval {$wig->mean} || undef;
}

sub draw_densefile {
  my $self = shift;
  my $feature = shift;
  my $densefile = shift;

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
  $self->create_parts_for_dense_feature($dense,$start,$end);
  $self->SUPER::draw(@_);
}

sub create_parts_for_dense_feature {
  my $self = shift;
  my ($dense,$start,$end) = @_;


#  my $span = $self->width;
  my $span = $self->scale> 1 ? $end - $start : $self->width;
  my $data = $dense->values($start,$end,$span);
  my $points_per_span = ($end-$start+1)/$span;
  my @parts;

  for (my $i=0; $i<$span;$i++) {
    my $offset = $i * $points_per_span;
    my $value  = shift @$data;
    next unless defined $value;
    push @parts,
      Bio::Graphics::Feature->new(-score => $value,
				  -start => int($start + $i * $points_per_span),
				  -end   => int($start + $i * $points_per_span));
  }
  $self->{parts} = [];
  $self->add_feature(@parts);
}

sub minmax {
  my $self  = shift;
  my $parts = shift;
  if (my $wig = $self->wig) {
    my $max = $self->option('max_score');
    my $min = $self->option('min_score');
    $max = $wig->max unless defined $max;
    $min = $wig->min unless defined $min;
    return ($min,$max);
  } else {
    return $self->SUPER::minmax($parts);
  }
}

sub subsample {
  my $self = shift;
  my ($data,$start,$span) = @_;
  my $points_per_span = @$data/$span;
  my @parts;
  for (my $i=0; $i<$span;$i++) {
    my $offset = $i * $points_per_span;
    my $value  = $data->[$offset + $points_per_span/2];
    push @parts,Bio::Graphics::Feature->new(-score => $value,
					    -start => int($start + $i * $points_per_span),
					    -end   => int($start + $i * $points_per_span));
  }
  return @parts;
}

sub create_parts_for_segment {
  my $self = shift;
  my ($seg,$start,$end) = @_;
  my $seg_start = $seg->start;
  my $seg_end   = $seg->end;
  my $step      = $seg->step;
  my $span      = $seg->span;

  # clip, because wig files do no clipping
  $seg_start = $start      if $seg_start < $start;
  $seg_end   = $end        if $seg_end   > $end;

  return unless $start < $end;

  # get data values across the area
  my @data = $seg->values($start,$end);

  # create a series of parts
  my @parts;
  for (my $i = $start; $i <= $end ; $i += $step) {
    my $data_point = shift @data;
    push @parts,Bio::Graphics::Feature->new(-score => $data_point,
					   -start => $i,
					   -end   => $i + $step - 1);
  }
  $self->{parts} = [];
  $self->add_feature(@parts);
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

Bio::Graphics::Glyph::wiggle_xyplot - An xyplot plot compatible with dense "wig"data

=head1 SYNOPSIS

  See <Bio::Graphics::Panel> and <Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph works like the regular xyplot but takes value data in
Bio::Graphics::Wiggle file format:

TODO! UPDATE DOCUMENTATION FOR DENSE FILES

 reference = chr1
 ChipCHIP Feature1 1..10000 wigfile=./test.wig
 ChipCHIP Feature2 10001..20000 wigfile=./test.wig
 ChipCHIP Feature3 25001..35000 wigfile=./test.wig

The "wigfile" attribute gives a relative or absolute pathname to a
Bio::Graphics::Wiggle format file. The data consist of a packed binary
representation of the values in the feature, using a constant step
such as present in tiling array data. Wigfiles are created using the
Bio::Graphics::Wiggle module or the wiggle2gff3.pl script, currently
both part of the gbrowse package.

=head2 OPTIONS

In addition to all the xyplot glyph options, the following options are
recognized:

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
