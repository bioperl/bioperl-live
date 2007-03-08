# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

use File::Spec;
use constant IMAGES => File::Spec->catfile(qw(t data biographics));
use constant FILES => File::Spec->catfile(qw(t data biographics));
use constant IMAGE_TESTS => 0;

my $error;

BEGIN { 
    $error = 0;
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 16 + (IMAGE_TESTS ? 3 : 0);
    plan tests => $NUMTESTS;
    
    require_ok('Bio::Graphics::FeatureFile');
    require_ok('Bio::Graphics');
}

exit 0 if $error;

SKIP: {
    eval {
        require GD;
        require Text::Shellwords;
    };
    if( $@ ) {
        skip("GD or Text::Shellwords modules are not installed. ".
        "This means that Bio::Graphics module is unusable. ".
        "Skipping tests.",$NUMTESTS-2);
    }
    my $verbose = -1;
    my $write   = 0;
    
    ## End of black magic.
    ##
    ## Insert additional test code below but remember to change
    ## the print "1..x\n" in the BEGIN block to reflect the
    ## total number of tests that will be run. 
    
    my @images = IMAGE_TESTS ? qw(t1 t2 t3) : ();
    
    # parse command line arguments
    while (@ARGV && $ARGV[0] =~ /^--?(\w+)/) {
      my $arg = $1;
      if ($arg eq 'write') {
        warn "Writing regression test images into ",IMAGES,".........\n";
        $write++;
      }
      shift;
    }
    
    
    foreach (@images) {
      if ($write) { warn "$_...\n"; do_write($_) } else { eval { do_compare($_) } }
    }
    
    my $data  = Bio::Graphics::FeatureFile->new(-file => FILES . "/feature_data.txt") or die;
    ok defined $data;
    is $data->render, 5;
    is $data->setting(general=>'pixels'), 750;
    is $data->setting('general'), 4;
    is $data->setting, 6;
    is $data->glyph('EST'), 'segments';
    
    my %style = $data->style('EST');
    is $style{-connector}, 'solid';
    is $style{-height}, 5;
    is $style{-bgcolor}, 'yellow';
    
    is $data->configured_types, 5;
    is @{$data->features('EST')}, 5;
    
    my $thing = $data->features('EST');
    
    my ($feature) = grep {$_->name eq 'Predicted gene 1'} @{$data->features('FGENESH')};
    ok $feature;
    is $feature->desc, "Pfam";
    is $feature->score, 20;
}

sub do_write {
  my $test = shift;
  my $canpng = GD::Image->can('png');
  my $output_file = IMAGES . ($canpng ? "/$test.png" : "/$test.gif");
  my $test_sub    = $test;
  my $panel       = eval "$test_sub()" or die "Couldn't run test: $@";
  open OUT,">$output_file" or die "Couldn't open $output_file for writing: $!";
  print OUT $canpng ? $panel->gd->png : $panel->gd->gif;
  close OUT;
}

sub do_compare {
  my $test = shift;
  my $canpng = GD::Image->can('png');
  my @input_files = glob(IMAGES . ($canpng ? "/$test/*.png" : "/$test/*.gif"));
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
  my $data  = Bio::Graphics::FeatureFile->new(-file => FILES . "/feature_data.txt") or die;
  my ($tracks,$panel) = $data->render;
  return $panel;
}
