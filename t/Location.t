# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 13 }
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::SimilarityPair;
use Bio::SeqFeature::FeaturePair;

ok(1);

my $simple = new Bio::Location::Simple('-start' => 10, '-end' => 20,
				       '-strand' => 1);
ok $simple->isa('Bio::LocationI') && $simple->isa('Bio::RangeI');

ok $simple->start, 10;
ok $simple->end, 20;

my $generic = new Bio::SeqFeature::Generic('-start' => 5, '-end' => 30, 
					   '-strand' => 1);

ok $generic->isa('Bio::SeqFeatureI') && $generic->isa('Bio::RangeI');
ok $generic->start, 5;
ok $generic->end, 30;

my $similarity = new Bio::SeqFeature::SimilarityPair();

my $feat1 = new Bio::SeqFeature::Generic('-start' => 30, '-end' => 43, 
					 '-strand' => -1);
my $feat2 = new Bio::SeqFeature::Generic('-start' => 80, '-end' => 90, 
					 '-strand' => -1);

my $featpair = new Bio::SeqFeature::FeaturePair('-feature1' => $feat1,
						'-feature2' => $feat2 );

my $feat3 = new Bio::SeqFeature::Generic('-start' => 35, '-end' => 50, 
					 '-strand' => -1);

ok $featpair->start, 30;
ok $featpair->end,  43;
ok $featpair->length, 14;

ok $featpair->overlaps($feat3);
ok $generic->overlaps($simple);
ok $generic->contains($simple);
