# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 45);
	
    use_ok('Bio::SeqFeature::Generic');
}

my @funcs = qw(start end length strand overlaps contains
    equals intersection union overlap_extent disconnected_ranges
    offsetStranded subtract);

my $i = 1;
while (my $func = shift @funcs ) {
    $i++;

    # test for presence of method
    ok exists $Bio::RangeI::{$func};
    
    # union get caught in an infinite loop w/o parameters; skip invoke test.
    next if $func eq 'union';
    
    # call to strand complains without a value; skip invoke test.
    next if $func eq 'disconnected_ranges';
    
    # test invocation of method
    eval { $Bio::RangeI::{$func}->(); };
    ok($@);
}

### unit tests for subtract method ###
# contributed by Stephen Montgomery (sm8 at sanger.ac.uk), who also
# wrote the subtract method
my $feature1 =  Bio::SeqFeature::Generic->new( -start => 1, -end =>
1000, -strand => 1);
my $feature2 =  Bio::SeqFeature::Generic->new( -start => 100, -end =>
900, -strand => -1);

my $subtracted = $feature1->subtract($feature2);
ok(defined($subtracted));
is(scalar(@$subtracted), 2);
foreach my $range (@$subtracted) {
    ok($range->start == 1 || $range->start == 901);
    ok($range->end == 99 || $range->end == 1000);
}

$subtracted = $feature2->subtract($feature1);
ok(!defined($subtracted));
$subtracted = $feature1->subtract($feature2, 'weak');
ok(!defined($subtracted));
$subtracted = $feature1->subtract($feature2, 'strong');
ok(!defined($subtracted));

my $feature3 =  Bio::SeqFeature::Generic->new( -start => 500, -end =>
1500, -strand => 1);
$subtracted = $feature1->subtract($feature3);
ok(defined($subtracted));
is scalar(@$subtracted), 1;
my $subtracted_i = @$subtracted[0];
is($subtracted_i->start, 1);
is($subtracted_i->end, 499);


# ---------------
# Added Bio::Location::SplitLocationI support to subtract().  --jhannah 20091018
$feature1 =  Bio::SeqFeature::Generic->new();
$feature2 =  Bio::SeqFeature::Generic->new();
my $loc = Bio::Location::Split->new();
$loc->add_sub_Location(Bio::Location::Simple->new(-start=>100, -end=>200, -strand=>1));
$loc->add_sub_Location(Bio::Location::Simple->new(-start=>300, -end=>400, -strand=>1));
$loc->add_sub_Location(Bio::Location::Simple->new(-start=>500, -end=>600, -strand=>1));
$feature1->location($loc);
$loc = Bio::Location::Split->new();
$loc->add_sub_Location(Bio::Location::Simple->new(-start=>350, -end=>400, -strand=>1));
$loc->add_sub_Location(Bio::Location::Simple->new(-start=>500, -end=>510, -strand=>1));
$feature2->location($loc);
$subtracted = $feature1->subtract($feature2);
is(@$subtracted, 3,                              "subtract() of split features");
is($subtracted->[0]->start, 100,                 "   0 start");
is($subtracted->[0]->end,   200,                 "   0 end");
is($subtracted->[1]->start, 300,                 "   1 start");
is($subtracted->[1]->end,   349,                 "   1 end");
is($subtracted->[2]->start, 511,                 "   2 start");
is($subtracted->[2]->end,   600,                 "   2 end");
# ---------------


