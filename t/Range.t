# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 16 }

use Bio::Range;

ok(1);

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

my ($start, $stop);

($start, $stop) = $range->union($range2);
ok ($start, 10);
ok ($stop, 25 );

($start, $stop) = $range->intersection($range2);
ok ($start, 15);
ok ($stop, 20);

ok !($range->contains($range2));
ok !($range2->contains($range));
ok ($range->overlaps($range2));
ok ($range2->overlaps($range));

ok ($range->strand(0), 0);
ok ($range->start(27), 27);
ok ($range->end(28), 28);

