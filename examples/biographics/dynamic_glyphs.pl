#!/usr/bin/perl

use lib '.','../..','./blib/lib','../../blib/lib','../..';
use strict;
use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

chomp (my $PKG = shift);
$PKG or die "\nUsage: lots_of_glyphs IMAGE_CLASS
\t- where IMAGE_CLASS is one of GD or GD::SVG
\t- GD generate png output; GD::SVG generates SVG.\n";

my $ftr = 'Bio::Graphics::Feature';

my $segment = $ftr->new(-start=>-100,-end=>1400,-name=>'ZK154',-type=>'clone');
my $zk154_1 = $ftr->new(-start=>-50,-end=>800,-name=>'ZK154.1',-type=>'gene');
my $zk154_2 = $ftr->new(-segments=>[[200,300],[380,800]],-name=>'ZK154.2',-type=>'gene');
my $zk154_3 = $ftr->new(-start=>900,-end=>1200,-name=>'ZK154.3',-type=>'gene');

my $zed_27 = $ftr->new(-segments=>[[550,600],[800,950],[1200,1300]],
		   -name=>'zed-27',
		   -subtype=>'exon',-type=>'transcript');
my $abc3 = $ftr->new(-segments=>[[100,200],[350,400],[500,550]],
		    -name=>'abc53',
		     -strand => -1,
		    -subtype=>'exon',-type=>'transcript');
my $xyz4 = $ftr->new(-segments=>[[40,80],[100,120],[200,280],[300,320]],
		     -name=>'xyz4',
		     -subtype=>'predicted',-type=>'alignment');

my $m3 = $ftr->new(-segments=>[[20,40],[30,60],[90,270],[290,300]],
		   -name=>'M3',
		   -subtype=>'predicted',-type=>'alignment');

my $bigone = $ftr->new(-segments=>[[-200,-120],[90,270],[290,300]],
		   -name=>'big one',
		   -subtype=>'predicted',-type=>'alignment');

my $fred_12 = $ftr->new(-segments=>[$xyz4,$zed_27],
			-type => 'group',
			-name =>'fred-12');

my $confirmed_exon1 = $ftr->new(-start=>1,-stop=>20,
				-type=>'exon',
				-source=>'confirmed',
				-name => 'confirmed1',
			       );
my $predicted_exon1 = $ftr->new(-start=>30,-stop=>50,
				-type=>'exon',
				-name=>'predicted1',
				-source=>'predicted');
my $predicted_exon2 = $ftr->new(-start=>60,-stop=>100,
				-name=>'predicted2',
				-type=>'exon',-source=>'predicted');

my $confirmed_exon3 = $ftr->new(-start=>150,-stop=>190,
				-type=>'exon',-source=>'confirmed',
			       -name=>'abc123');
my $partial_gene = $ftr->new(-segments=>[$confirmed_exon1,$predicted_exon1,$predicted_exon2,$confirmed_exon3],
			     -name => 'partial gene',
			     -type => 'transcript',
			     -source => '(from a big annotation pipeline)'
			    );
my @segments = $partial_gene->segments;
my $score = 10;
foreach (@segments) {
  $_->score($score);
  $score += 10;
}

my $panel = Bio::Graphics::Panel->new(
				      -gridcolor => 'lightcyan',
				      -grid => 1,
				      -segment => $segment,
				      -spacing => 15,
				      -width   => 600,
				      -pad_top  => 20,
				      -pad_bottom  => 20,
				      -pad_left => 20,
				      -pad_right=> 20,
				      -key_style => 'between',
         		 	      -image_class=> $PKG,
				     );
my @colors = $panel->color_names();

my $t = $panel->add_track(
			  transcript2 => [$abc3,$zed_27],
			  -label => 1,
			  -bump => 1,
			  -key => 'Prophecies',
			  #		  -tkcolor => $colors[rand @colors],
			 );
$t->configure(-bump=>1);
$panel->add_track($segment,
		  -glyph => 'arrow',
		  -label => sub {scalar localtime},
#		  -labelfont => 'gdMediumBoldFont',
		  -double => 1,
		  -bump => 0,
		  -height => 10,
		  -arrowstyle=>'regular',
		  -linewidth=>1,
		  -tick => 2,
		 );

$panel->add_track(generic => [$segment,$abc3,$zk154_1,[$zk154_2,$xyz4]],
		  -label     => sub { $_[-1]->level == 0 } ,
		  -bgcolor   => sub { shift->primary_tag eq 'predicted' ? 'green' : 'blue'},
		  -connector => sub { my $primary_tag = shift->primary_tag;
				      $primary_tag eq 'transcript' ? 'hat'
				    : $primary_tag eq 'alignment'  ? 'solid'
				    : 'solid'},
		  -connector_color => 'black',
		  -height => 10,
		  -bump => 1,
#		  -tkcolor => $colors[rand @colors],
		  -key => 'Signals',
		 );

my $track = $panel->add_track('transcript2'=> [$bigone],
			      -label   => 1,
			      -connector => 'solid',
			      -point  => 0,
			      -orient => 'N',
			      -height => 8,
			      -base => 1,
			      -relative_coords => 1,
			      -tick  => 2,
			      -bgcolor => 'red',
			      -key     => 'Dynamically Added');
#$track->add_feature($bigone,$zed_27,$abc3);
#$track->add_group($predicted_exon1,$predicted_exon2,$confirmed_exon3);
$track->add_group($bigone,$zed_27,$zk154_2,$bigone);

my $gd    = $panel->gd;
my @boxes = $panel->boxes;
my $red   = $panel->translate_color('red');
for my $box (@boxes) {
  my ($feature,@points) = @$box;
}
my $type = ($PKG eq 'GD') ? 'png' : 'svg';
print $gd->$type;
