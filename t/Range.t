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
    plan tests => 17;
}

use Bio::Range;

ok 1;

my $range = Bio::Range->new(-start=>10,
                            -end=>20,
			    -strand=>1);
ok defined $range;
ok $range->strand, 1;

my $range2 = Bio::Range->new(-start=>15,
                             -end=>25,
			     -strand=>1);

ok defined $range2;
ok $range2->strand, 1;

my $r = Bio::Range->new();
ok $r->strand(0), 0;
ok $r->start(27), 27;
ok $r->end(28), 28;

ok $r->intersection($range2), undef;

$r = $range->union($range2);
ok $r->start, 10;
ok $r->end, 25;

$r = $range->intersection($range2);
ok $r->start, 15;
ok $r->end, 20;

ok !($range->contains($range2));
ok !($range2->contains($range));
ok $range->overlaps($range2);
ok $range2->overlaps($range);


