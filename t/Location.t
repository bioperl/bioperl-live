# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    eval { require Test::More; };
    if( $@ ) {
		use lib 't/lib';
    }
    use Test::More;
    plan tests => 104;
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Location::Split');
    use_ok('Bio::Location::Fuzzy');
    
    use_ok('Bio::SeqFeature::Generic');
    use_ok('Bio::SeqFeature::SimilarityPair');
    use_ok('Bio::SeqFeature::FeaturePair');
}

my $simple = new Bio::Location::Simple('-start' => 10, '-end' => 20,
				       '-strand' => 1, -seq_id => 'my1');
isa_ok($simple, 'Bio::LocationI');
isa_ok($simple, 'Bio::RangeI');

is($simple->start, 10, 'Bio::Location::Simple tests');
is($simple->end, 20);
is($simple->seq_id, 'my1');

my ($loc) = $simple->each_Location();
ok($loc);
is("$loc", "$simple");

my $generic = new Bio::SeqFeature::Generic('-start' => 5, '-end' => 30, 
					   '-strand' => 1);

isa_ok($generic,'Bio::SeqFeatureI', 'Bio::SeqFeature::Generic' );
isa_ok($generic,'Bio::RangeI');
is($generic->start, 5);
is($generic->end, 30);

my $similarity = new Bio::SeqFeature::SimilarityPair();

my $feat1 = new Bio::SeqFeature::Generic('-start' => 30, '-end' => 43, 
					 '-strand' => -1);
my $feat2 = new Bio::SeqFeature::Generic('-start' => 80, '-end' => 90, 
					 '-strand' => -1);

my $featpair = new Bio::SeqFeature::FeaturePair('-feature1' => $feat1,
						'-feature2' => $feat2 );

my $feat3 = new Bio::SeqFeature::Generic('-start' => 35, '-end' => 50, 
					 '-strand' => -1);

is($featpair->start, 30,'Bio::SeqFeature::FeaturePair tests');
is($featpair->end,  43);

is($featpair->length, 14);

ok($featpair->overlaps($feat3));
ok($generic->overlaps($simple), 'Bio::SeqFeature::Generic tests');
ok($generic->contains($simple));

# fuzzy location tests
my $fuzzy = new Bio::Location::Fuzzy('-start'  =>'<10', 
				     '-end'    => 20,
				     -strand   =>1, 
				     -seq_id   =>'my2');

is($fuzzy->strand, 1, 'Bio::Location::Fuzzy tests');
is($fuzzy->start, 10);
is($fuzzy->end,20);
ok(!defined $fuzzy->min_start);
is($fuzzy->max_start, 10);
is($fuzzy->min_end, 20);
is($fuzzy->max_end, 20);
is($fuzzy->location_type, 'EXACT');
is($fuzzy->start_pos_type, 'BEFORE');
is($fuzzy->end_pos_type, 'EXACT');
is($fuzzy->seq_id, 'my2');
is($fuzzy->seq_id('my3'), 'my3');

($loc) = $fuzzy->each_Location();
ok($loc);
is("$loc", "$fuzzy");

# split location tests
my $splitlocation = new Bio::Location::Split;
my $f = new Bio::Location::Simple(-start  => 13,
				  -end    => 30,
				  -strand => 1);
$splitlocation->add_sub_Location($f);
is($f->start, 13, 'Bio::Location::Split tests');
is($f->min_start, 13);
is($f->max_start,13);


$f = new Bio::Location::Simple(-start  =>30,
			       -end    =>90,
			       -strand =>1);
$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Simple(-start  =>18,
			       -end    =>22,
			       -strand =>1);
$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Simple(-start  =>19,
			       -end    =>20,
			       -strand =>1);

$splitlocation->add_sub_Location($f);

$f = new Bio::Location::Fuzzy(-start  =>"<50",
			      -end    =>61,
			      -strand =>1);
is($f->start, 50);
ok(! defined $f->min_start);
is($f->max_start, 50);

is(scalar($splitlocation->each_Location()), 4);

$splitlocation->add_sub_Location($f);

is($splitlocation->max_end, 90);
is($splitlocation->min_start, 13);
is($splitlocation->end, 90);
is($splitlocation->start, 13);
is($splitlocation->sub_Location(),5);


is($fuzzy->to_FTstring(), '<10..20');
$fuzzy->strand(-1);
is($fuzzy->to_FTstring(), 'complement(<10..20)');
is($simple->to_FTstring(), '10..20');
$simple->strand(-1);
is($simple->to_FTstring(), 'complement(10..20)');
is( $splitlocation->to_FTstring(), 
    'join(13..30,30..90,18..22,19..20,<50..61)');

# test for bug #1074
$f = new Bio::Location::Simple(-start  => 5,
			       -end    => 12,
			       -strand => -1);
$splitlocation->add_sub_Location($f);
is( $splitlocation->to_FTstring(), 
    'join(13..30,30..90,18..22,19..20,<50..61,complement(5..12))',
	'Bugfix 1074');
$splitlocation->strand(-1);
is( $splitlocation->to_FTstring(), 
    'complement(join(13..30,30..90,18..22,19..20,<50..61,5..12))');

$f = new Bio::Location::Fuzzy(-start => '45.60',
			      -end   => '75^80');

is($f->to_FTstring(), '(45.60)..(75^80)');
$f->start('20>');
is($f->to_FTstring(), '>20..(75^80)');

# test that even when end < start that length is always positive

$f = new Bio::Location::Simple(-verbose => -1,
			       -start   => 100, 
			       -end     => 20, 
			       -strand  => 1);

is($f->length, 81, 'Positive length');
is($f->strand,-1);

# test that can call seq_id() on a split location;
$splitlocation = new Bio::Location::Split(-seq_id => 'mysplit1');
is($splitlocation->seq_id,'mysplit1', 'seq_id() on Bio::Location::Split');
is($splitlocation->seq_id('mysplit2'),'mysplit2');


# Test Bio::Location::Exact

ok(my $exact = new Bio::Location::Simple(-start    => 10, 
					 -end      => 20,
					 -strand   => 1, 
					 -seq_id   => 'my1'));
isa_ok($exact, 'Bio::LocationI');
isa_ok($exact, 'Bio::RangeI');

is( $exact->start, 10, 'Bio::Location::Simple EXACT');
is( $exact->end, 20);
is( $exact->seq_id, 'my1');
is( $exact->length, 11);
is( $exact->location_type, 'EXACT');

ok ($exact = new Bio::Location::Simple(-start         => 10, 
				      -end           => 11,
				      -location_type => 'IN-BETWEEN',
				      -strand        => 1, 
				      -seq_id        => 'my2'));

is($exact->start, 10, 'Bio::Location::Simple IN-BETWEEN');
is($exact->end, 11);
is($exact->seq_id, 'my2');
is($exact->length, 0);
is($exact->location_type, 'IN-BETWEEN');

eval {
    $exact = new Bio::Location::Simple(-start         => 10, 
				       -end           => 12,
				       -location_type => 'IN-BETWEEN');
};
ok( $@, 'Testing error handling' );

# testing error when assigning 10^11 simple location into fuzzy
eval {
    ok $fuzzy = new Bio::Location::Fuzzy(-start         => 10, 
					 -end           => 11,
					 -location_type => '^',
					 -strand        => 1, 
					 -seq_id        => 'my2');
};
ok( $@ );

$fuzzy = new Bio::Location::Fuzzy(-location_type => '^',
				  -strand        => 1, 
				  -seq_id        => 'my2');

$fuzzy->start(10);
eval { $fuzzy->end(11) };
ok($@);

$fuzzy = new Bio::Location::Fuzzy(-location_type => '^',
				  -strand        => 1, 
				  -seq_id        =>'my2');

$fuzzy->end(11);
eval {
    $fuzzy->start(10);
};
ok($@);

# testing coodinate policy modules

use_ok('Bio::Location::WidestCoordPolicy');
use_ok('Bio::Location::NarrowestCoordPolicy');
use_ok('Bio::Location::AvWithinCoordPolicy');

$f = new Bio::Location::Fuzzy(-start => '40.60',
			      -end   => '80.100');
is $f->start, 40, 'Default coodinate policy';
is $f->end, 100;
is $f->length, 61;
is $f->to_FTstring, '(40.60)..(80.100)';
isa_ok($f->coordinate_policy, 'Bio::Location::WidestCoordPolicy');

# this gives an odd location string; is it legal?
$f->coordinate_policy(new Bio::Location::NarrowestCoordPolicy);
is $f->start, 60, 'Narrowest coodinate policy';
is $f->end, 80;
is $f->length, 21;
is $f->to_FTstring, '(60.60)..(80.80)';
isa_ok($f->coordinate_policy, 'Bio::Location::NarrowestCoordPolicy');

# this gives an odd location string
$f->coordinate_policy(new Bio::Location::AvWithinCoordPolicy);
is $f->start, 50, 'Average coodinate policy';
is $f->end, 90;
is $f->length, 41;
is $f->to_FTstring, '(50.60)..(80.90)';
isa_ok($f->coordinate_policy, 'Bio::Location::AvWithinCoordPolicy');

# to complete the circle
$f->coordinate_policy(new Bio::Location::WidestCoordPolicy);
is $f->start, 40, 'Widest coodinate policy';
is $f->end, 100;
is $f->length, 61;
is $f->to_FTstring, '(40.60)..(80.100)';
isa_ok($f->coordinate_policy, 'Bio::Location::WidestCoordPolicy');

