package Bio::Graphics::Glyph::ternary_plot;
# Draw ternary (triangle plots)

use strict;
use base qw(Bio::Graphics::Glyph::generic);
use Bio::Graphics::Glyph::xyplot;

use constant Sin60 =>0.866025403784439;
use constant Tan60 =>1.73205080756888;

sub calculate_side {
  my $self = shift;
  $self->option('height') / Sin60;
}

# positioning this properly is a bit tricky, because if the side of the triangle
# is greater than the width of the feature, then we need to add extra left and right
# padding.
sub pad_left {
  my $self  = shift;
  my $feature = $self->feature;
  my $left   = $self->SUPER::pad_left;

  my $side   = $self->calculate_side;
  my ($a,$b) = $self->map_pt($feature->start,$feature->stop);
  my $width  = abs($b-$a);
  my $extra  = $width > $side ? 0 : ($side-$width)/2;

  return $extra > $left ? $extra : $left;
}

sub pad_right {
  my $self = shift;
  my $right = $self->SUPER::pad_right;
  my $side   = $self->calculate_side;
  my $feature = $self->feature;
  my ($a,$b) = $self->map_pt($feature->start,$feature->stop);
  my $width  = abs($b-$a);
  my $extra  = $width > $side ? 0 : ($side-$width)/2;
  return $extra > $right ? $extra : $right;
}

sub pad_top {
  my $self = shift;
  my $pad  = $self->SUPER::pad_top;
  return $pad unless $self->option('vertices');
  my $font = $self->image_class->gdTinyFont();
  my $lh = $font->height/2;
  return $pad > $lh ? $pad : $lh;
}

sub pad_bottom {
  my $self = shift;
  my $pad  = $self->SUPER::pad_bottom;
  return $pad unless $self->option('vertices');
  my $font = $self->image_class->gdTinyFont();
  my $lh = $font->height/2;
  return $pad > $lh ? $pad : $lh;
}

sub triples {
  my $self = shift;
  my $triples = $self->option('triples');
  return $triples if defined $triples;
  my @triples = $self->feature->get_tag_values('triples');
  for my $t (@triples) {
    next if ref $t && ref $t eq 'ARRAY';  # already in right format
    $t = [split /[,\w]/,$t];
  }
  return \@triples;
}

sub draw_component {
  my $self = shift;
  my $gd = shift;
  my $fg = $self->fgcolor;
  my $bg = $self->bgcolor;

  # position the left (A) edge
  my ($x1,$y1,$x2,$y2) = $self->calculate_boundaries(@_);
  my $xmid = ($x1+$x2)/2;
  my $side = $self->calculate_side;

  my $left  = $xmid - $side/2;
  my $right = $left + $side;
  my $top   = $y1;
  my $bottom=$y2;

  # draw the triangle
  my $poly_pkg = $self->polygon_package;
  my $poly     = $poly_pkg->new();
  $poly->addPt($left,$bottom);
  $poly->addPt($right,$bottom);
  $poly->addPt($xmid,$top);
  $gd->polygon($poly,$fg);

  # draw vertex labels, if any
  my $fontcolor  = $self->fontcolor;
  my $font = $self->image_class->gdTinyFont();
  my $lh = $font->height;
  my $lw = $font->width;

  if (my $vertex_labels = $self->option('vertices')) {
    my @labels = @$vertex_labels;
    $gd->string($font,$left-$lw*length($labels[0]),$bottom+$lh/2,$labels[0],$fontcolor);
    $gd->string($font,$right,$bottom+$lh/2,$labels[1],$fontcolor);
    $gd->string($font,$xmid-$lw*length($labels[2])-3,$top-3,$labels[2],$fontcolor);
  }

  # get triples
  my $data = $self->triples;
  for my $triple (@$data) {
    my ($a,$b,$c,$color,$label) = @$triple;
    $color = defined $color ? $self->factory->translate_color($color) : $bg;
    my $x = $xmid    + $side * ($b - $a)/2;
    my $y = $bottom  - $side * ($c * Tan60)/2;
    draw_disc($gd,$x,$y,3,$color);
    if ($label) {
      $gd->string($font,$x+3,$y,$label,$fontcolor);
    }
  }
}

sub draw_disc {
  my ($gd,$x,$y,$pr,$color) = @_;
  $gd->filledArc($x,$y,$pr,$pr,0,360,$color);
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::ternary_plot - Draw ternary plot data

=head1 SYNOPSIS

 #!/usr/bin/perl

 use strict;
 use warnings;

 use Bio::Graphics;
 use Bio::Graphics::Feature;

 my $segment  = Bio::Graphics::Feature->new(-start=>1,-end=>700);
 my $snp1     = Bio::Graphics::Feature->new(-start     => 500,
 					   -end       => 501,
					   -name      => 'rs000001',
					   -attributes=> {triples => [
								      [0.01, 0.81, 0.18, 'red',  'CEPH'],
								      [0.25, 0.25, 0.50, 'blue', 'JPT+CHB'],
								      [0.81, 0.01, 0.18, 'green','YRI'],
								     ]
							  }
					  );
 my $snp2     = Bio::Graphics::Feature->new(-start     => 300,
					   -end       => 301,
					   -name      => 'rs12345',
					   -attributes=> {triples => [
								      [0.04, 0.64, 0.32, 'red',  'Controls'],
								      [0.16, 0.36, 0.48, 'blue', 'Cases'],
								     ]
							  }
					  );

 my $panel = Bio::Graphics::Panel->new(-segment=>$segment,-width=>800);

 $panel->add_track($segment,-glyph=>'arrow',-double=>1,-tick=>2);
 $panel->add_track([$snp1,$snp2],
		  -glyph    => 'ternary_plot',
		  -height   => 80,
		  -fgcolor  => 'lightgrey',
		  -vertices => ['AA','GG','AG'],
		  -label    => 1,
		 );

 print $panel->png;

=head1 DESCRIPTION

This glyph draws a light gray equilateral triangle with its base
centered on the feature. The top of the equilateral triangle is equal
to the specified height. To look good, please choose a height of E<gt>=
15.

Inside, the glyph will plot one or more data points using ternary plot
conventions (see http://en.wikipedia.org/wiki/Ternary_plot). The data
consists of a series of (A,B,C) triplets chosen such that the range of
each component is [0.0,1.0] and A + B + C = 1.0. The left, right and
apex of the triangle represent the proportions of A, B and C
respectively. As a component approaches 1.0, it gets closer to its
corresponding vertex.

The data can be represented as one or more feature tags called "triples"
each in the format:

   A1,B1,C1,<color>,<label>   # (color and label are optional)

or as a callback specified by the option B<-triples>, which should
return a list of arrays, where each array is a triple, followed by an
optional color. E.G.

 sub my_calback {
   my $feature = shift;
   return [[0.1,0.5,0.4,'red','pt1'],[0.2,0.2,0.6,'blue','pt2'],[0.8,0.2,0.0,'green','pt4']];
 }

The color, if it is missing, will be the same as the bgcolor.

=head2 OPTIONS

In addition to the common options, the following glyph-specific
options are recognized:

  Option      Description 
  ------      -----------

  -triples    The callback to return triple data.
  -vertices   Labels for the left,right & top vertices

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
L<Bio::Graphics::Glyph::whiskerplot>,
L<Bio::DB::GFF>,
L<Bio::SeqI>,
L<Bio::SeqFeatureI>,
L<Bio::Das>,
L<GD>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut
