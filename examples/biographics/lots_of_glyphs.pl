#!/usr/bin/perl -w

use lib '.','..','./blib/lib','../blib/lib';
use strict;

use Bio::Graphics::Panel;
use Bio::Graphics::Feature;

my $ftr = 'Bio::Graphics::Feature';

my $segment = $ftr->new(-start=>-100,-end=>1000,-name=>'ZK154',-type=>'clone');
my $zk154_1 = $ftr->new(-start=>-50,-end=>800,-name=>'ZK154.1',-type=>'gene');
my $zk154_2 = $ftr->new(-start=>380,-end=>500,-name=>'ZK154.2',-type=>'gene');
my $zk154_3 = $ftr->new(-start=>900,-end=>1200,-name=>'ZK154.3',-type=>'gene');

my $zed_27 = $ftr->new(-segments=>[[400,500],[550,600],[800,950]],
		   -name=>'zed-27',
		   -subtype=>'exon',-type=>'gene');
my $abc3 = $ftr->new(-segments=>[[100,200],[350,400],[500,550]],
		    -name=>'abc53',
		     -strand => -1,
		    -subtype=>'exon',-type=>'gene');
my $xyz4 = $ftr->new(-segments=>[[40,80],[100,120],[200,280],[300,320]],
		     -name=>'xyz4',
		     -subtype=>'predicted',-type=>'alignment');

my $m3 = $ftr->new(-segments=>[[20,40],[30,60],[90,270],[290,300]],
		   -name=>'M3',
		   -subtype=>'predicted',-type=>'alignment');

my $bigone = $ftr->new(-segments=>[[-200,-120],[90,270],[290,300]],
		       -name=>'big one',
		       -subtype=>'predicted',-type=>'gene');

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
#				      -grid => [50,100,150,200,250,300,310,320,330],
				      -gridcolor => 'lightcyan',
				      -grid => 1,
				      -segment => $segment,
#				      -offset => 300,
#				      -length  => 1000,
				      -spacing => 15,
				      -width   => 600,
				      -pad_top  => 20,
				      -pad_bottom  => 20,
				      -pad_left => 20,
				      -pad_right=> 20,
#				      -bgcolor => 'teal',
#				      -key_style => 'between',
				      -key_style => 'bottom',
				     );
my @colors = $panel->color_names();

my $t = $panel->add_track(
			  #		  generic => [$abc3,$zed_27],
			  transcript2 => [$abc3,$zed_27],
			  -label => 1,
			  -bump => 1,
			  -key => 'Prophecies',
			  #		  -tkcolor => $colors[rand @colors],
			 );
$t->configure(-bump=>1);
$panel->add_track($segment,
		  -glyph => 'arrow',
		  -label => 'base pairs',
		  -double => 1,
		  -bump => 0,
		  -height => 10,
		  -arrowstyle=>'regular',
		  -linewidth=>1,
#		  -tkcolor => $colors[rand @colors],
		  -tick => 2,
		 );
$panel->unshift_track(generic => [$segment,$zk154_1,$zk154_2,$zk154_3,[$xyz4,$zed_27]],
		      -label     => sub { my $feature = shift; $feature->sub_SeqFeature>0},
		      -bgcolor   => sub { shift->primary_tag eq 'predicted' ? 'olive' : 'red'},
		      -connector => sub { my $feature = shift;
					  my $type = $feature->primary_tag;
					  $type eq 'group'      ? 'dashed'
					    : $type eq 'transcript' ? 'hat'
					      : $type eq 'alignment'  ? 'solid'
						: undef},
		      -all_callbacks => 1,
		      -connector_color => 'black',
		      -height => 10,
		      -bump => 1,
		      -linewidth=>2,
		      #		  -tkcolor => $colors[rand @colors],
		      -key => 'Signs',
		 );

my $track = $panel->add_track(-glyph=> sub { shift->primary_tag eq 'gene' ? 'transcript2': 'generic'},
			      -label   => sub { $_[-1]->level == 0 } ,
			      -connector => sub { return shift->type eq 'group' ? 'dashed' : ''},
			      -point  => 0,
			      -orient => 'N',
			      -height => 8,
			      -base => 1,
			      -relative_coords => 1,
			      -tick  => 2,
			      -all_callbacks => 1,
			      -bgcolor => 'red',
			      -key     => 'Dynamically Added');
$track->add_feature($bigone,$zed_27,$abc3);
$track->add_group($predicted_exon1,$predicted_exon2,$confirmed_exon3);

$panel->add_track(
		  [$abc3,$zed_27,$partial_gene],
		  -bgcolor   => sub { shift->source_tag eq 'predicted' ? 'green' : 'blue'},
		  -glyph   => 'transcript',
#		  -glyph   => sub { my $feature = shift; 
#				    return $feature->source_tag eq 'predicted'
#				      ? 'ellipse' : 'transcript'},
		  -label       => sub { shift->sub_SeqFeature > 0 },
#		  -label       => 1,
#		  -description => sub { shift->sub_SeqFeature > 0 },
		  -description => sub {
		    my $feature = shift;
		    return 1   if $feature->primary_tag eq 'transcript';
		    return '*' if $feature->source_tag eq 'predicted';
		    return;
		  },
		  -font2color  => 'red',
		  -bump => +1,
#		  -tkcolor => $colors[rand @colors],
		  -key => 'Portents',
		 );
$panel->add_track(segments => [$segment,$zk154_1,[$zk154_2,$xyz4]],
		  -label     => 1,
		  -bgcolor   => sub { shift->primary_tag eq 'predicted' ? 'green' : 'blue'},
		  -connector => sub { my $primary_tag = shift->primary_tag;
				      $primary_tag eq 'transcript' ? 'hat'
				    : $primary_tag eq 'alignment'  ? 'solid'
				    : undef},
		  -connector_color => 'black',
		  -height => 10,
		  -bump => 1,
#		  -tkcolor => $colors[rand @colors],
		  -key => 'Signals',
		 );
$panel->add_track(generic => [],
		  -key => 'Foobar');

$panel->add_track(graded_segments => $partial_gene,
		  -bgcolor =>'blue',
		  -label   => 1,
		  -key     => 'Scored thing');

$panel->add_track(diamond => [$segment,$zk154_1,$zk154_2,$zk154_3,$xyz4,$zed_27],
		  -bgcolor =>'blue',
		  -label   => 1,
		  -key     => 'pointy thing');

#print $panel->png;

my $gd    = $panel->gd;
my @boxes = $panel->boxes;
my $red   = $panel->translate_color('red');
for my $box (@boxes) {
  my ($feature,@points) = @$box;
#  $gd->rectangle(@points,$red);
}
#$gd->filledRectangle(0,0,20,200,1);
#$gd->filledRectangle(600-20,0,600,200,1);
print $gd->png;

