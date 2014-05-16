# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 109);
	
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Location::Split');
    use_ok('Bio::Location::Fuzzy');
    use_ok('Bio::SeqFeature::Generic');
    use_ok('Bio::SeqFeature::SimilarityPair');
    use_ok('Bio::SeqFeature::FeaturePair');
    use_ok('Bio::SeqFeature::Lite');
}

my $simple = Bio::Location::Simple->new('-start' => 10, '-end' => 20,
				       '-strand' => 1, -seq_id => 'my1');
isa_ok($simple, 'Bio::LocationI');
isa_ok($simple, 'Bio::RangeI');

is($simple->start, 10, 'Bio::Location::Simple tests');
is($simple->end, 20);
is($simple->seq_id, 'my1');

my ($loc) = $simple->each_Location();
ok($loc);
is("$loc", "$simple");

my $generic = Bio::SeqFeature::Generic->new('-start' => 5, '-end' => 30, 
					   '-strand' => 1);

isa_ok($generic,'Bio::SeqFeatureI', 'Bio::SeqFeature::Generic' );
isa_ok($generic,'Bio::RangeI');
is($generic->start, 5);
is($generic->end, 30);

my $lite_pos  = Bio::SeqFeature::Lite->new(-start => 1000, -stop => 2000, -strand => '+');
my $lite_neg  = Bio::SeqFeature::Lite->new(-start => 1000, -end  => 2000, -strand => '-');
my $lite_none = Bio::SeqFeature::Lite->new(-start => 1000, -stop => 2000, -strand => '.');
is($lite_pos->strand,  1);
is($lite_neg->strand, -1);
is($lite_neg->end,  2000);
is($lite_neg->stop, 2000);
is($lite_none->strand, 0);

my $similarity = Bio::SeqFeature::SimilarityPair->new();

my $feat1 = Bio::SeqFeature::Generic->new('-start' => 30, '-end' => 43, 
					 '-strand' => -1);
my $feat2 = Bio::SeqFeature::Generic->new('-start' => 80, '-end' => 90, 
					 '-strand' => -1);

my $featpair = Bio::SeqFeature::FeaturePair->new('-feature1' => $feat1,
						'-feature2' => $feat2 );

my $feat3 = Bio::SeqFeature::Generic->new('-start' => 35, '-end' => 50, 
					 '-strand' => -1);

is($featpair->start, 30,'Bio::SeqFeature::FeaturePair tests');
is($featpair->end,  43);

is($featpair->length, 14);

ok($featpair->overlaps($feat3));
ok($generic->overlaps($simple), 'Bio::SeqFeature::Generic tests');
ok($generic->contains($simple));

# fuzzy location tests
my $fuzzy = Bio::Location::Fuzzy->new('-start'  =>'<10', 
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
my $splitlocation = Bio::Location::Split->new();
my $f = Bio::Location::Simple->new(-start  => 13,
				  -end    => 30,
				  -strand => 1);
$splitlocation->add_sub_Location($f);
is($f->start, 13, 'Bio::Location::Split tests');
is($f->min_start, 13);
is($f->max_start,13);


$f = Bio::Location::Simple->new(-start  =>30,
			       -end    =>90,
			       -strand =>1);
$splitlocation->add_sub_Location($f);

$f = Bio::Location::Simple->new(-start  =>18,
			       -end    =>22,
			       -strand =>1);
$splitlocation->add_sub_Location($f);

$f = Bio::Location::Simple->new(-start  =>19,
			       -end    =>20,
			       -strand =>1);

$splitlocation->add_sub_Location($f);

$f = Bio::Location::Fuzzy->new(-start  =>"<50",
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
$f = Bio::Location::Simple->new(-start  => 5,
			       -end    => 12,
			       -strand => -1);
$splitlocation->add_sub_Location($f);
is( $splitlocation->to_FTstring(), 
    'join(13..30,30..90,18..22,19..20,<50..61,complement(5..12))',
	'Bugfix 1074');
$splitlocation->strand(-1);
is( $splitlocation->to_FTstring(), 
    'complement(join(13..30,30..90,18..22,19..20,<50..61,5..12))');

$f = Bio::Location::Fuzzy->new(-start => '45.60',
			      -end   => '75^80');

is($f->to_FTstring(), '(45.60)..(75^80)');
$f->start('20>');
is($f->to_FTstring(), '>20..(75^80)');

# test that even when end < start that length is always positive

$f = Bio::Location::Simple->new(-verbose => -1,
			       -start   => 100, 
			       -end     => 20, 
			       -strand  => 1);

is($f->length, 81, 'Positive length');
is($f->strand,-1);

# test that can call seq_id() on a split location;
$splitlocation = Bio::Location::Split->new(-seq_id => 'mysplit1');
is($splitlocation->seq_id,'mysplit1', 'seq_id() on Bio::Location::Split');
is($splitlocation->seq_id('mysplit2'),'mysplit2');


# Test Bio::Location::Exact

ok(my $exact = Bio::Location::Simple->new(-start    => 10, 
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

ok ($exact = Bio::Location::Simple->new(-start         => 10, 
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
    $exact = Bio::Location::Simple->new(-start         => 10, 
				       -end           => 12,
				       -location_type => 'IN-BETWEEN');
};
ok( $@, 'Testing error handling' );

# testing error when assigning 10^11 simple location into fuzzy
eval {
    ok $fuzzy = Bio::Location::Fuzzy->new(-start         => 10, 
					 -end           => 11,
					 -location_type => '^',
					 -strand        => 1, 
					 -seq_id        => 'my2');
};
ok( $@ );

$fuzzy = Bio::Location::Fuzzy->new(-location_type => '^',
				  -strand        => 1, 
				  -seq_id        => 'my2');

$fuzzy->start(10);
eval { $fuzzy->end(11) };
ok($@);

$fuzzy = Bio::Location::Fuzzy->new(-location_type => '^',
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

$f = Bio::Location::Fuzzy->new(-start => '40.60',
			      -end   => '80.100');
is $f->start, 40, 'Default coodinate policy';
is $f->end, 100;
is $f->length, 61;
is $f->to_FTstring, '(40.60)..(80.100)';
isa_ok($f->coordinate_policy, 'Bio::Location::WidestCoordPolicy');

# this gives an odd location string; is it legal?
$f->coordinate_policy(Bio::Location::NarrowestCoordPolicy->new());
is $f->start, 60, 'Narrowest coodinate policy';
is $f->end, 80;
is $f->length, 21;
is $f->to_FTstring, '(60.60)..(80.80)';
isa_ok($f->coordinate_policy, 'Bio::Location::NarrowestCoordPolicy');

# this gives an odd location string
$f->coordinate_policy(Bio::Location::AvWithinCoordPolicy->new());
is $f->start, 50, 'Average coodinate policy';
is $f->end, 90;
is $f->length, 41;
is $f->to_FTstring, '(50.60)..(80.90)';
isa_ok($f->coordinate_policy, 'Bio::Location::AvWithinCoordPolicy');

# to complete the circle
$f->coordinate_policy(Bio::Location::WidestCoordPolicy->new());
is $f->start, 40, 'Widest coodinate policy';
is $f->end, 100;
is $f->length, 61;
is $f->to_FTstring, '(40.60)..(80.100)';
isa_ok($f->coordinate_policy, 'Bio::Location::WidestCoordPolicy');
