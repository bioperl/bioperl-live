#!/usr/bin/perl

# file: render_blast1.pl
# This is code example 1 in the Graphics-HOWTO
# Author: Lincoln Stein

use strict;
use lib "$ENV{HOME}/projects/bioperl-live";
use Bio::Graphics;
use Bio::SeqFeature::Generic;

my $panel = Bio::Graphics::Panel->new(-length => 1000,-width  => 800);
my $track = $panel->add_track(-glyph => 'generic',-label  => 1);

while (<>) { # read blast file
  chomp;
  next if /^\#/;  # ignore comments
  my($name,$score,$start,$end) = split /\t+/;
  my $feature = Bio::SeqFeature::Generic->new(-seq_id=>$name,-score=>$score,
					      -start=>$start,-end=>$end);
  $track->add_feature($feature);
}

print $panel->png;

