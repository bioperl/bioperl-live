package Bio::Graphics::Glyph::whiskerplot;

use strict;
use base qw(Bio::Graphics::Glyph::xyplot);

sub lookup_draw_method {
  my $self = shift;
  return 'draw_whiskerplot';
}

sub draw_whiskerplot {
  my $self = shift;
  my ($gd,$left,$top) = @_;
  my @parts   = $self->parts;
  my $fgcolor = $self->fgcolor;

  for my $part (@parts) {
    my ($x1,$y1,$x2,$y2) = $part->calculate_boundaries($left,$top);

    my $y = $part->{_y_position};

    # get the range and tendency
    if (my $range =  $part->get_range) {

      unshift @$range,undef if @$range == 4;  # backward compatibility

      my ($median,$range_low,$range_high,$lower_quartile,$higher_quartile) = @$range;
      $y = $part->{_y_position} = $self->score2position($median) if defined $median;

      # draw the quartile box
      my ($box_top,$box_bottom) = ($y,$y);
      if (defined $lower_quartile && defined $higher_quartile) {
	$box_top    = $self->score2position($higher_quartile);
	$box_bottom = $self->score2position($lower_quartile);
	$self->filled_box($gd,$x1,$box_top,$x2,$box_bottom);
      }

      # calculate positions of the range whiskers
      if (defined $range_low && defined $range_high) {
	my $range_top    = $self->score2position($range_high);
	my $range_bottom = $self->score2position($range_low);
	my $center       = ($x1+$x2)/2;
	my $whisker_left  = $center-5;
	my $whisker_right = $center+5;
	$whisker_left     = $x1 if $whisker_left  < $x1;
	$whisker_right    = $x2 if $whisker_right > $x2;

	# top whisker
	$gd->line($center,$box_top,$center,$range_top,$fgcolor);
	$gd->line($whisker_left,$range_top,$whisker_right,$range_top,$fgcolor);

	# bottom whisker
	$gd->line($center,$box_bottom,$center,$range_bottom,$fgcolor);
	$gd->line($whisker_left,$range_bottom,$whisker_right,$range_bottom,$fgcolor);
      }
    }

    # draw the median
    $gd->line($x1,$y,$x2,$y,$fgcolor);

  }
}

sub get_range {
  my $self  = shift;
  my $range = $self->option('range');
  return $range if defined $range;
  # otherwise get it from the feature
  return [$self->feature->get_tag_values('range')];
}

1;

__END__

=head1 NAME

Bio::Graphics::Glyph::whiskerplot - The whiskerplot glyph

=head1 SYNOPSIS

  See L<Bio::Graphics::Panel> and L<Bio::Graphics::Glyph>.

=head1 DESCRIPTION

This glyph is used for drawing features associated with numeric data
using "box and whisker" style data points, which display the mean
value, extreme ranges and first and third quartiles (or standard
deviation). The boxes drawn by this glyph are similar to
L<http://www.abs.gov.au/websitedbs/D3310116.NSF/0/3c35ac1e828c23ef4a2567ac0020ec8a?OpenDocument>,
except that they are oriented vertically so that the position and
height of the box indicates the mean value and spread of the data, and
the width indicates the genomic extent of the value.

Like the xyplot glyph (from which it inherits the whiskerplot is
designed to work on a single feature group that contains subfeatures.
It is the subfeatures that carry the score information. The best way
to arrange for this is to create an aggregator for the feature.  We'll
take as an example a histogram of repeat density in which interval are
spaced every megabase and the score indicates the number of repeats in
the interval; we'll assume that the database has been loaded in in
such a way that each interval is a distinct feature with the method
name "density" and the source name "repeat".  Furthermore, all the
repeat features are grouped together into a single group (the name of
the group is irrelevant).  If you are using Bio::DB::GFF and
Bio::Graphics directly, the sequence of events would look like this:

  my $agg = Bio::DB::GFF::Aggregator->new(-method    => 'repeat_density',
                                          -sub_parts => 'density:repeat');
  my $db  = Bio::DB::GFF->new(-dsn=>'my_database',
                              -aggregators => $agg);
  my $segment  = $db->segment('Chr1');
  my @features = $segment->features('repeat_density');

  my $panel = Bio::Graphics::Panel->new;
  $panel->add_track(\@features,
                    -glyph => 'xyplot',
                    -scale => 'both',
);

If you are using Generic Genome Browser, you will add this to the
configuration file:

  aggregators = repeat_density{density:repeat}
                clone alignment etc

Note that it is a good idea to add some padding to the left and right
of the panel; otherwise the scale will be partially cut off by the
edge of the image.

The mean (or median) of the data will be taken from the feature
score. The range and quartile data must either be provided in a
feature tag named "range", or must be generated dynamically by a
-range callback option passed to add_track. The data returned by the
tag or option should be an array reference containing the following
five fields:

 [$median,$range_low,$range_high,$quartile_low,$quartile_high]

where $range_low and $range_high correspond to the low and high value
of the "whiskers" and $quartile_low and $quartile_high correspond to
the low and high value of the "box."

If $median is undef or missing, then the score field of the feature
will be used instead. It may be useful to repeat the median in the
score field in any case, in order to allow the minimum and maximum
range calculations of the graph itself to occur.

See Examples for three ways of generating an image.


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

  -range        Callback to return median,    none - data comes from feature "range" tag
                range and quartiles for each
                sub feature

=head1 EXAMPLES

Here are three examples of how to use this glyph.

=head2 Example 1: Incorporating the numeric data in each subfeature

 #!/usr/bin/perl
 use strict;

 use Bio::Graphics;
 use Bio::SeqFeature::Generic;

 my $bsg = 'Bio::SeqFeature::Generic';

 my $feature = $bsg->new(-start=>0,-end=>1000);

 for (my $i=0;$i<1000;$i+=20) {
   my $y = (($i-500)/10)**2;
   my $range = make_range($y);
   my $part = $bsg->new(-start=>$i,-end=>$i+16,
 		       -score=>$y,-tag => { range=>$range });
   $feature->add_SeqFeature($part);
 }

 my $panel = Bio::Graphics::Panel->new(-length=>1000,-width=>800,-key_style=>'between',
 				      -pad_left=>40,-pad_right=>40);
 $panel->add_track($feature,
 		  -glyph=>'arrow',
 		  -double=>1,
 		  -tick=>2);

 $panel->add_track($feature,
 		  -glyph=>'whiskerplot',
 		  -scale=>'both',
 		  -height=>200,
 		  -min_score => -500,
 		  -key  =>'Whiskers',
 		  -bgcolor => 'orange',
 		 );
 print $panel->png;

 sub make_range {
   my $score        = shift;
   my $range_top    = $score + 5*sqrt($score) + rand(50);
   my $range_bottom = $score - 5*sqrt($score) - rand(50);
   my $quartile_top    = $score + 2*sqrt($score) + rand(50);
   my $quartile_bottom = $score - 2*sqrt($score) - rand(50);
   return [$score,$range_bottom,$range_top,$quartile_bottom,$quartile_top];
 }

=head2 Example 2: Generating the range data with a callback

 #!/usr/bin/perl
 use strict;

 use Bio::Graphics;
 use Bio::SeqFeature::Generic;

 my $bsg = 'Bio::SeqFeature::Generic';
 my $feature = $bsg->new(-start=>0,-end=>1000);

 for (my $i=0;$i<1000;$i+=20) {
   my $y = (($i-500)/10)**2;
   my $part = $bsg->new(-start=>$i,-end=>$i+16,-score=>$y);
   $feature->add_SeqFeature($part);
 }

 my $panel = Bio::Graphics::Panel->new(-length=>1000,-width=>800,-key_style=>'between',
 				      -pad_left=>40,-pad_right=>40);
 $panel->add_track($feature,
 		  -glyph=>'arrow',
 		  -double=>1,
 		  -tick=>2);

 $panel->add_track($feature,
 		  -glyph=>'whiskerplot',
 		  -scale=>'both',
 		  -height=>200,
 		  -min_score => -500,
 		  -key  =>'Whiskers',
 		  -bgcolor => 'orange',
 		  -range => \&make_range,
 		 );
 print $panel->png;

 sub make_range {
   my $feature = shift;
   my $score        = $feature->score;
   my $range_top    = $score + 5*sqrt($score) + rand(50);
   my $range_bottom = $score - 5*sqrt($score) - rand(50);
   my $quartile_top    = $score + 2*sqrt($score) + rand(50);
   my $quartile_bottom = $score - 2*sqrt($score) - rand(50);
   return [$score,$range_bottom,$range_top,$quartile_bottom,$quartile_top];
 }

=head2 Example 3: Generating the image from a FeatureFile

=over 4

=item The file:

 [general]
 pixels = 840
 pad_left = 40
 pad_right = 40

 [contig]
 glyph     = arrow
 double    = 1
 tick      = 2

 [data]
 glyph     = whiskerplot
 scale     = both
 height    = 200
 min_score = -500
 max_score = 2800
 key       = Whiskers
 bgcolor   = orange

 chr1	.	contig	1	1000	.	.	.	Contig chr1
 chr1	.	data	0	16	2500	.	.	Dataset data1; range 2209,2769,2368,2619
 chr1	.	data	20	36	2304	.	.	Dataset data1; range 2051,2553,2163,2435
 chr1	.	data	40	56	2116	.	.	Dataset data1; range 1861,2384,1983,2253
 chr1	.	data	60	76	1936	.	.	Dataset data1; range 1706,2181,1819,2059
 chr1	.	data	80	96	1764	.	.	Dataset data1; range 1516,1995,1646,1849
 chr1	.	data	100	116	1600	.	.	Dataset data1; range 1359,1834,1513,1699
 chr1	.	data	120	136	1444	.	.	Dataset data1; range 1228,1654,1330,1565
 chr1	.	data	140	156	1296	.	.	Dataset data1; range 1105,1520,1198,1385
 chr1	.	data	160	176	1156	.	.	Dataset data1; range 983,1373,1062,1270
 chr1	.	data	180	196	1024	.	.	Dataset data1; range 853,1184,914,1116
 chr1	.	data	200	216	900	.	.	Dataset data1; range 722,1093,801,965
 chr1	.	data	220	236	784	.	.	Dataset data1; range 621,945,724,859
 chr1	.	data	240	256	676	.	.	Dataset data1; range 532,833,605,742
 chr1	.	data	260	276	576	.	.	Dataset data1; range 433,714,485,653
 chr1	.	data	280	296	484	.	.	Dataset data1; range 331,600,418,545
 chr1	.	data	300	316	400	.	.	Dataset data1; range 275,535,336,459
 chr1	.	data	320	336	324	.	.	Dataset data1; range 198,434,270,374
 chr1	.	data	340	356	256	.	.	Dataset data1; range 167,378,219,322
 chr1	.	data	360	376	196	.	.	Dataset data1; range 114,303,118,249
 chr1	.	data	380	396	144	.	.	Dataset data1; range 39,248,87,197
 chr1	.	data	400	416	100	.	.	Dataset data1; range 17,173,68,141
 chr1	.	data	420	436	64	.	.	Dataset data1; range -14,125,18,84
 chr1	.	data	440	456	36	.	.	Dataset data1; range -8,74,11,64
 chr1	.	data	460	476	16	.	.	Dataset data1; range -46,77,0,43
 chr1	.	data	480	496	4	.	.	Dataset data1; range -40,43,-7,36
 chr1	.	data	500	516	0	.	.	Dataset data1; range -43,0,-43,22
 chr1	.	data	520	536	4	.	.	Dataset data1; range -6,52,-4,54
 chr1	.	data	540	556	16	.	.	Dataset data1; range -5,38,-27,52
 chr1	.	data	560	576	36	.	.	Dataset data1; range -43,109,18,66
 chr1	.	data	580	596	64	.	.	Dataset data1; range -1,134,3,112
 chr1	.	data	600	616	100	.	.	Dataset data1; range 49,186,69,124
 chr1	.	data	620	636	144	.	.	Dataset data1; range 79,225,71,169
 chr1	.	data	640	656	196	.	.	Dataset data1; range 124,289,120,266
 chr1	.	data	660	676	256	.	.	Dataset data1; range 154,378,197,320
 chr1	.	data	680	696	324	.	.	Dataset data1; range 220,439,249,396
 chr1	.	data	700	716	400	.	.	Dataset data1; range 291,511,331,458
 chr1	.	data	720	736	484	.	.	Dataset data1; range 350,627,400,572
 chr1	.	data	740	756	576	.	.	Dataset data1; range 446,718,502,633
 chr1	.	data	760	776	676	.	.	Dataset data1; range 515,833,576,777
 chr1	.	data	780	796	784	.	.	Dataset data1; range 606,959,724,856
 chr1	.	data	800	816	900	.	.	Dataset data1; range 747,1058,799,1004
 chr1	.	data	820	836	1024	.	.	Dataset data1; range 817,1231,958,1089
 chr1	.	data	840	856	1156	.	.	Dataset data1; range 961,1341,1069,1225
 chr1	.	data	860	876	1296	.	.	Dataset data1; range 1103,1511,1219,1385
 chr1	.	data	880	896	1444	.	.	Dataset data1; range 1218,1660,1338,1535
 chr1	.	data	900	916	1600	.	.	Dataset data1; range 1377,1828,1496,1703
 chr1	.	data	920	936	1764	.	.	Dataset data1; range 1547,2020,1674,1858
 chr1	.	data	940	956	1936	.	.	Dataset data1; range 1691,2188,1824,2043
 chr1	.	data	960	976	2116	.	.	Dataset data1; range 1869,2376,2019,2225
 chr1	.	data	980	996	2304	.	.	Dataset data1; range 2040,2554,2178,2418

=item The script to render it

 #!/usr/bin/perl

 use strict;
 use Bio::Graphics::FeatureFile;

 my $data = Bio::Graphics::FeatureFile->new(-file=>'test.gff');

 my(undef,$panel) = $data->render;
 print $panel->png;

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

Copyright (c) 2001 Cold Spring Harbor Laboratory

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

