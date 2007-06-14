# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw(@funcs);

BEGIN {

    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if ($@) {
        use lib 't';
    }
    use Test;
    @funcs = qw(start end length strand overlaps contains
        equals intersection union overlap_extent disconnected_ranges
        offsetStranded subtract);
    plan tests => 37;
}

use Bio::RangeI;

my $i = 1;
my $func;
while ( $func = shift @funcs ) {
    $i++;

    # test for presence of method
    if ( exists $Bio::RangeI::{$func} ) {
        ok(1);


      # union get caught in an infinite loop w/o parameters; skip invoke test.
        next if $func eq 'union';

        # call to strand complains without a value; skip invoke test.
        next if $func eq 'disconnected_ranges';

        # test invocation of method
        eval { $Bio::RangeI::{$func}->(); };
        ok($@);
    }
    else {
        ok(0);
    }
}

### unit tests for subtract method ###
# contributed by Stephen Montgomery (sm8 at sanger.ac.uk), who also
# wrote the subtract method
use Bio::SeqFeature::Generic;
use Data::Dumper;

my $feature1 =  Bio::SeqFeature::Generic->new( -start => 1, -end =>
1000, -strand => 1);
my $feature2 =  Bio::SeqFeature::Generic->new( -start => 100, -end =>
900, -strand => -1);

my $subtracted = $feature1->subtract($feature2);
ok(defined($subtracted));
ok(scalar(@$subtracted) == 2);
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
ok(scalar(@$subtracted) == 1);
my $subtracted_i = @$subtracted[0];
ok($subtracted_i->start == 1);
ok($subtracted_i->end == 499);

# end