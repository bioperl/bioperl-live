#!/usr/bin/perl

# file: render_blast4.pl
# This is code example 4 in the Graphics-HOWTO
# Author: Lincoln Stein

use strict;
use lib "$ENV{HOME}/projects/bioperl-live";
use Bio::Graphics;
use Bio::SearchIO;

my $file = shift or die "Usage: render4.pl <blast file>\n";

my $searchio = Bio::SearchIO->new(-file   => $file,
				  -format => 'blast') or die "parse failed";


my $result = $searchio->next_result() or die "no result";

my $panel = Bio::Graphics::Panel->new(-length    => $result->query_length,
				      -width     => 800,
				      -pad_left  => 10,
				      -pad_right => 10,
				     );

my $full_length = Bio::SeqFeature::Generic->new(-start=>1,-end=>$result->query_length,
						-seq_id=>$result->query_name);
$panel->add_track($full_length,
		  -glyph   => 'arrow',
		  -tick    => 2,
		  -fgcolor => 'black',
		  -double  => 1,
		  -label   => 1,
		 );

my $track = $panel->add_track(-glyph       => 'graded_segments',
			      -label       => 1,
			      -connector   => 'dashed',
			      -bgcolor     => 'blue',
			      -font2color  => 'red',
			      -sort_order  => 'high_score',
			      -description => sub {
				my $feature = shift;
				return unless $feature->has_tag('description');
				my ($description) = $feature->each_tag_value('description');
				my $score = $feature->score;
				"$description, score=$score";
			       });

while( my $hit = $result->next_hit ) {
  next unless $hit->significance < 1E-20;
  my $feature = Bio::SeqFeature::Generic->new(-score   => $hit->raw_score,
					      -seq_id   => $hit->name,
					      -tag     => {
							   description => $hit->description
							  },
					     );
  while( my $hsp = $hit->next_hsp ) {
    $feature->add_sub_SeqFeature($hsp,'EXPAND');
  }

  $track->add_feature($feature);
}

print $panel->png;

