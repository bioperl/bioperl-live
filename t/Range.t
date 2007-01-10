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
	plan tests => 44;
}

use Bio::Range;

ok(1);

my $range = Bio::Range->new(-start=>10,
                            -end=>20,
									 -strand=>1);
ok(defined $range);
ok($range->strand, 1);

my $range2 = Bio::Range->new(-start=>15,
                             -end=>25,
									  -strand=>1);

ok(defined $range2);
ok($range2->strand, 1);

my $r = Bio::Range->new();
ok($r->strand(0), 0);
ok($r->start(27), 27);
ok($r->end(28), 28);

ok(! defined $r->intersection($range2));

$r = $range->union($range2);
ok($r->start, 10);
ok($r->end, 25);

$r = $range->intersection($range2);
ok($r->start, 15);
ok($r->end, 20);
ok($r->strand, 1);

# intersection and union can also take lists
my $range3 = Bio::Range->new(-start=>18,-end=>30);
$r = $range->intersection([$range2, $range3]);
ok($r->start == 18 && $r->end == 20);
$r = Bio::Range->intersection([$range, $range2, $range3]);
ok($r->start == 18 && $r->end == 20);
$r = $range->union($range2, $range3);
ok($r->start == 10 && $r->end == 30);
$r = Bio::Range->union($range, $range2, $range3);
ok($r->start == 10 && $r->end == 30);
$range3->start(21);
ok(! defined $range->intersection([$range2, $range3]));

ok !($range->contains($range2));
ok !($range2->contains($range));
ok ($range->overlaps($range2));
ok ($range2->overlaps($range));

# testing strand
$range3 = Bio::Range->new(-start => 15,
                             -end => 25,
									  -strand => 1);

my $range4 = Bio::Range->new(-start => 15,
									  -end => 25,
									  -strand => -1);

my $range5 = Bio::Range->new(-start => 15,
                             -end => 25,
									  -strand => 0);

my $range6 = Bio::Range->new(-start => 20,
									  -end => 30,
									  -strand => -1);

ok $range3->_ignore($range4);     # 1 & -1
ok $range3->_weak($range3);       # 1 & 1 true
ok $range3->_weak($range5);       # 1 & 0 true
ok ! ($range3->_weak($range4));   # 1 & -1 false
ok $range3->_strong($range3);     # 1 & 1 true
ok ! ($range3->_strong($range5)); # 1 & 0 false
ok ! ($range3->_strong($range4)); # 1 & -1 false

ok ! ( $range3->overlaps($range4,'weak'));
ok ! ( $range4->overlaps($range3,'weak'));
ok ! ( $range3->overlaps($range4,'strong')); 
ok ! ( $range4->overlaps($range3,'strong')); 

$range3->strand(0);

ok  ( $range3->overlaps($range4,'weak'));
ok  ( $range4->overlaps($range3,'weak')); 
ok ! ( $range3->overlaps($range4,'strong'));
ok ! ( $range4->overlaps($range3,'strong')); 

# if strands are different then intersection() should return 0...
$r = $range3->intersection($range4);
ok($r->strand, 0);

# or if both strands are -1 then -1 should be returned
$r = $range6->intersection($range4);
ok($r->strand, -1);

#test implemention of offsetStranded:
$r = Bio::Range->new(-start => 30, -end => 40, -strand => -1);
ok ($r->offsetStranded(-5,10)->toString, '(20, 45) strand=-1');
ok ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=-1');
$r->strand(1);
ok ($r->offsetStranded(-5,10)->toString, '(25, 50) strand=1');
ok ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=1');
