# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use File::Spec;

# In order to properly run the image comparison tests the images may need to be
# regenerated from scratch; this is primarily due to changes in GD versions, OS,
# Bio::Graphics, problems with storing binary data in CVS, etc.

# We'll need to reconfigure these tests to allow do_write()
# the ability to regenerate those files when passing the --write option

# for now, the image tests are turned off

use constant IMAGE_TESTS => 0;

BEGIN { 
  use lib '.';
  use Bio::Root::Test;
  
  test_begin(-tests => 48 + (IMAGE_TESTS ? 3 : 0),
             -requires_modules => [qw(GD)]);
  
  use_ok('Bio::Graphics::FeatureFile');
  use_ok('Bio::Graphics::Panel');
}

my $images = test_input_file('biographics');

my @images = IMAGE_TESTS ? qw(t1 t2 t3) : ();

# parse command line arguments
my $write = 0;
while (@ARGV && $ARGV[0] =~ /^--?(\w+)/) {
  my $arg = $1;
  if ($arg eq 'write') {
    warn "Writing regression test images into ",$images,".........\n";
    $write++;
  }
  shift;
}

foreach (@images) {
  if ($write) { warn "$_...\n"; do_write($_) } else { eval { do_compare($_) } }
}

my $data  = Bio::Graphics::FeatureFile->new(-file => test_input_file('biographics', 'feature_data.txt'),
					    -safe => 0,
    ) or die;
ok defined $data;
is $data->render, 5;
is $data->setting(general=>'pixels'), 750;
is $data->setting('general'), 3;
is $data->setting, 6;
is $data->glyph('EST'), 'segments';

my %style = $data->style('EST');
is $style{-connector}, 'solid';
is $style{-height}, 5;
is $style{-bgcolor}, 'yellow';

is $data->configured_types, 5;
is @{$data->features('EST')}, 5;

my $thing = $data->features('EST');
is $thing->[0]->seq_id,'B0511';

my ($feature) = grep {$_->name eq 'Predicted gene 1'} @{$data->features('FGENESH')};
ok $feature;
is $feature->desc, "Pfam";
is $feature->score, 20;

# test handling of things that look like comments
is $data->setting(EST=>'bgcolor'),'yellow';
is $data->setting(EST=>'fgcolor'),'#EE00FF';
is $data->setting(EST=>'link'),'http://www.google.com/search?q=$name#results';

# test handling of adding features
$data->add_type(TEST=>{bgcolor=>'green',
		       feature=>'test_feature',
		       glyph => 'generic'});
is $data->setting(TEST=>'bgcolor'),'green';
is $data->setting(TEST=>'feature'),'test_feature';
$data->add_feature(Bio::Graphics::FeatureBase->new(-seq_id    => 'chr1',
						   -start     => 1,
						   -end       => 1000,
						   -primary_tag=> 'test_feature'));
$data->add_feature(Bio::Graphics::FeatureBase->new(-seq_id    => 'chr2',
						   -start     => 2,
						   -end       => 2000,
						   -primary_tag=> 'test_feature'));
$data->add_feature(Bio::Graphics::FeatureBase->new(-seq_id    => 'chr3',
						   -start     => 3,
						   -end       => 3000),
		   'test_feature');
my @f = $data->features('test_feature');
is scalar @f,3;

# test FeatureBase
my $bfg   = 'Bio::Graphics::FeatureBase';
$feature  = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>1);
is $feature->seq_id,'chr2';
is $feature->start,201;
is $feature->end,300;
is $feature->strand,1;

# plus strand feature, plus strand ref sequence
my $ref   = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>1);
$feature->refseq($ref);
is $feature->start,1;
is $feature->end,100;
is $feature->strand,1;
is $feature->abs_start,201;
is $feature->abs_end,300;
is $feature->abs_strand,1;

# plus strand feature, minus strand ref sequence
$ref      = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>-1);
$feature->refseq($ref);
is $feature->start,100;   # expect flipping so that start > end
is $feature->end,1;
is $feature->strand,-1;

# minus strand feature, plus strand ref
$feature = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>-1);
$ref   = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>1);
$feature->refseq($ref);
is $feature->start,1;
is $feature->end,100;
is $feature->strand,-1;

# minus strand feature, minus strand ref
$ref      = $bfg->new(-seq_id=>'chr2',-start=>201,-end=>300,-strand=>-1);
$feature->refseq($ref);
is $feature->start,100;   # expect flipping so that start > end
is $feature->end,1;
is $feature->strand,1;

# test safety of callbacks
is $data->safe,0;
is ref $data->setting(SwissProt=>'fill'),'';
is eval{ref $data->code_setting(SwissProt=>'fill')},undef;

$data  = Bio::Graphics::FeatureFile->new(-file => test_input_file('biographics', 'feature_data.txt'),
					 -safe => 1,
    ) or die;
is $data->safe,1;
is ref $data->setting(SwissProt=>'fill'),'CODE';
is eval{ref $data->code_setting(SwissProt=>'fill')},'CODE';

exit 0;

sub do_write {
  my $test = shift;
  my $canpng = GD::Image->can('png');
  my $cangif = GD::Image->can('gif');
  my $test_sub    = $test;
  if ($canpng) {
    my $output_file = test_input_files($images,$test,"png");
    my $panel       = eval "$test_sub()" or die "Couldn't run test: $@";
    open OUT,">$output_file" or die "Couldn't open $output_file for writing: $!";
    print OUT $panel->gd->png;
    close OUT;
  }
  if ($cangif) {
    my $output_file = test_input_files($images, $test,"gif");
    my $panel       = eval "$test_sub()" or die "Couldn't run test: $@";
    open OUT,">$output_file" or die "Couldn't open $output_file for writing: $!";
    print OUT $panel->gd->gif;
    close OUT;
  }
}

sub test_input_files {
    my ($dir,$base,$suffix) = @_;
    my $index = 1;
    my $file;
    while (1) {
	$file = File::Spec->catfile($dir,$base,"version${index}.${suffix}");
	last unless -e $file;
	$index++;
    }
    return $file;
}

sub do_compare {
  my $test = shift;
  my $canpng = GD::Image->can('png');
  my @input_files = glob($images . ($canpng ? "/$test/*.png" : "/$test/*.gif"));
  my $test_sub    = $test;
  my $panel       = eval "$test_sub()" or die "Couldn't run test";
  my $ok = 0;
  my $test_data = $canpng ? $panel->gd->png : $panel->gd->gif;
  foreach (@input_files) {
    my $reference_data = read_file($_);
    if ($reference_data eq $test_data) {
      $ok++;
      last;
    }
  }
  ok($ok);
}

sub read_file {
  my $f = shift;
  open F,$f or die "Can't open $f: $!";
  binmode(F);
  my $data = '';
  while (read(F,$data,1024,length $data)) { 1 }
  close F;
  $data;
}


sub t1 {

  my $ftr = 'Bio::Graphics::Feature';

  my $segment = $ftr->new(-start=>1,-end=>1000,-name=>'ZK154',-type=>'clone');
  my $subseg1 = $ftr->new(-start=>1,-end=>500,-name=>'seg1',-type=>'gene');
  my $subseg2 = $ftr->new(-start=>250,-end=>500,-name=>'seg2',-type=>'gene');
  my $subseg3 = $ftr->new(-start=>250,-end=>500,-name=>'seg3',-type=>'gene');
  my $subseg4 = $ftr->new(-start=>1,-end=>400,-name=>'seg4',-type=>'gene');
  my $subseg5 = $ftr->new(-start=>400,-end=>800,-name=>'seg5',-type=>'gene');
  my $subseg6 = $ftr->new(-start=>550,-end=>800,-name=>'seg6',-type=>'gene');
  my $subseg7 = $ftr->new(-start=>550,-end=>800,-name=>'seg7',-type=>'gene');
  my $subseg8 = $ftr->new(-segments=>[[100,200],[300,400],[420,800]],-name=>'seg8',-type=>'gene');

  my $panel = Bio::Graphics::Panel->new(
					-grid => 1,
					-segment => $segment,
					-key_style => 'bottom');
  $panel->add_track(segments=>[$subseg1,$subseg2,$subseg3,$subseg4,
			       $subseg5,$subseg6,$subseg7,$subseg8],
		    -bump => 1,
		    -label => 1,
		    -key => '+1 bumping');
  $panel->add_track(segments=>[$subseg1,$subseg2,$subseg3,$subseg4,
			       $subseg5,$subseg6,$subseg7,$subseg8],
		    -bump => -1,
		    -label => 1,
		    -bgcolor => 'blue',
		    -key => '-1 bumping');
  $panel->add_track(segments=>[$subseg1,$subseg2,$subseg3,$subseg4,
			       $subseg5,$subseg6,$subseg7,$subseg8],
		    -bump => +2,
		    -label => 1,
		    -bgcolor => 'orange',
		    -key => '+2 bumping');
  $panel->add_track(segments=>[$subseg1,$subseg2,$subseg3,$subseg4,
			       $subseg5,$subseg6,$subseg7,$subseg8],
		    -bump => -2,
		    -label => 1,
		    -bgcolor => 'yellow',
		    -key => '-2 bumping');
  return $panel;
}


sub t2 {
  my $ftr = 'Bio::Graphics::Feature';

  my $segment = $ftr->new(-start=>-100,-end=>1000,-name=>'ZK154',-type=>'clone');
  my $zk154_1 = $ftr->new(-start=>-50,-end=>800,-name=>'ZK154.1',-type=>'gene');
  my $zk154_2 = $ftr->new(-start=>380,-end=>500,-name=>'ZK154.2',-type=>'gene');
  my $zk154_3 = $ftr->new(-start=>900,-end=>1200,-name=>'ZK154.3',-type=>'gene');

  my $zed_27 = $ftr->new(-segments=>[[400,500],[550,600],[800,950]],
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
				  -desc=>'confirmed',
				  -name => 'confirmed1',
				 );
  my $predicted_exon1 = $ftr->new(-start=>30,-stop=>50,
				  -type=>'exon',
				  -name=>'predicted1',
				  -desc=>'predicted');
  my $predicted_exon2 = $ftr->new(-start=>60,-stop=>100,
				  -name=>'predicted2',
				  -type=>'exon',-desc=>'predicted');

  my $confirmed_exon3 = $ftr->new(-start=>150,-stop=>190,
				  -type=>'exon',-desc=>'confirmed',
				  -name=>'abc123');
  my $partial_gene = $ftr->new(-segments=>[$confirmed_exon1,$predicted_exon1,$predicted_exon2,$confirmed_exon3],
			       -name => 'partial gene',
			       -type => 'transcript',
			       -desc => '(from a big annotation pipeline)'
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
					-empty_tracks => 'suppress',
				       );
  my @colors = $panel->color_names();

  my $t = $panel->add_track(
			    transcript2 => [$abc3,$zed_27],
			    -label => 1,
			    -bump => 1,
			    -key => 'Prophecies',
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
			-key => 'Signs',
			-empty_tracks => 'suppress',
		       );

  my $track = $panel->add_track(-glyph=> sub { shift->primary_tag =~ /transcript|alignment/ ? 'transcript2': 'generic'},
				-label   => sub { $_[-1]->level == 0 } ,
				-connector => sub { return shift->type eq 'group' ? 'dashed' : 'hat'},
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
		    -label       => sub { shift->sub_SeqFeature > 0 },
		    -description => sub {
		      my $feature = shift;
		      return 1   if $feature->primary_tag eq 'transcript';
		      return '*' if $feature->source_tag eq 'predicted';
		      return;
		    },
		    -font2color  => 'red',
		    -bump => +1,
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
		    -key => 'Signals',
		   );
  $panel->add_track(generic => [],
		    -key => 'Empty');

  $panel->add_track(graded_segments => $partial_gene,
		    -bgcolor =>'blue',
		    -vary_fg => 1,
		    -label   => 1,
		    -key     => 'Scored thing');

  $panel->add_track(diamond => [$segment,$zk154_1,$zk154_2,$zk154_3,$xyz4,$zed_27],
		    -bgcolor =>'blue',
		    -label   => 1,
		    -key     => 'pointy thing');
  return $panel;
}

sub t3 {
  my $data  = Bio::Graphics::FeatureFile->new(-file => test_input_file('biographics', 'feature_data.txt')) or die;
  my ($tracks,$panel) = $data->render;
  return $panel;
}
