# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 32; 
}

use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Location::Fuzzy;

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

my $splitlocation = new Bio::Location::Split;
$splitlocation->add_sub_Location(new Bio::Location::Simple('-start'=>1,
							   '-end'=>30,
							   '-strand'=>1));
$splitlocation->add_sub_Location(new Bio::Location::Simple('-start'=>50,
							   '-end'=>61,
							   '-strand'=>1));

ok($splitlocation->max_end, 61);
ok($splitlocation->min_start, 1);

ok($splitlocation->sub_Location(),2);

my $fuzzy = new Bio::Location::Fuzzy('-start'=>'<10', '-end' => 20, 
				     -strand=>1);
				     
ok($fuzzy->start, 10);
ok($fuzzy->end,20);
ok($fuzzy->fuzzy_start, '<10');
ok($fuzzy->fuzzy_end, 20);
ok($fuzzy->fuzzy_range, '..');

my ($encode,$pt) = $fuzzy->_fuzzypoint('>5');
ok($encode, -2);
ok($pt, 5);

($encode,$pt) = $fuzzy->_fuzzypoint('<5');
ok($encode, -1);
ok($pt, 5);

($encode,$pt) = $fuzzy->_fuzzypoint('5>');
ok($encode, 1);
ok($pt, 5);

($encode,$pt) = $fuzzy->_fuzzypoint('5<');
ok($encode, 2);
ok($pt, 5);

($encode,$pt) = $fuzzy->_fuzzypoint('5');
ok($encode, 0);
ok($pt, 5);

$fuzzy->verbose(-1);
($encode,$pt) = $fuzzy->_fuzzypoint('badstr');
ok($encode, undef);
ok($pt, undef);
$fuzzy->verbose(0);
