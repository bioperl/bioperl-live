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
    plan tests => 43;
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

# fuzzy location tests
my $fuzzy = new Bio::Location::Fuzzy('-start' =>'<10', '-end' => 20, 
				     -strand=>1);
				     
ok($fuzzy->strand, 1);
ok($fuzzy->start, 10);
ok($fuzzy->end,20);
ok($fuzzy->min_start, undef);
ok($fuzzy->max_start, 10);
ok($fuzzy->min_end, 20);
ok($fuzzy->max_end, 20);
ok($fuzzy->loc_type, 'EXACT');
ok($fuzzy->start_pos_type, 'BEFORE');
ok($fuzzy->end_pos_type, 'EXACT');

# split location tests
my $splitlocation = new Bio::Location::Split;
my $f = new Bio::Location::Simple('-start'=>13,
				  '-end'=>30,
				  '-strand'=>1);
$splitlocation->add_sub_Location($f);
ok($f->start, 13);
ok($f->min_start, 13);
ok($f->max_start,13);

$f = new Bio::Location::Simple('-start'=>30,
			       '-end'=>90,
			       '-strand'=>1);
$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Simple('-start'=>18,
			       '-end'=>22,
			       '-strand'=>1);
$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Simple('-start'=>19,
				  '-end'=>20,
			       '-strand'=>1);

$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Fuzzy('-start'=>"<50",
			      '-end'=>61,
			      '-strand'=>1);
ok($f->start, 50);
ok($f->min_start, undef);
ok($f->max_start, 50);

$splitlocation->add_sub_Location($f);

ok($splitlocation->max_end, 90);
ok($splitlocation->min_start, 13);
ok($splitlocation->end, 90);
ok($splitlocation->start, 13);
ok($splitlocation->sub_Location(),5);

ok($fuzzy->to_FTstring(), '<10..20');
$fuzzy->strand(-1);
ok($fuzzy->to_FTstring(), 'complement(<10..20)');
ok($simple->to_FTstring(), '10..20');
$simple->strand(-1);
ok($simple->to_FTstring(), 'complement(10..20)');
ok( $splitlocation->to_FTstring(), 'join(13..30,30..90,18..22,19..20,<50..61)');
$f = new Bio::Location::Fuzzy(-start => '45.60',
			      -end   => '75^80');

ok($f->to_FTstring(), '45.60..75^80');
$f->start('20>');
ok($f->to_FTstring(), '>20..75^80');

# test that even when end < start that length is always positive

$f = new Bio::Location::Simple(-verbose => -1,-start => 100, -end => 20, -strand => 1);

ok($f->length, 81);
ok($f->strand, -1);
