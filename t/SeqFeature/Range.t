# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 49);
	
    use_ok('Bio::Range');
}

my $range = Bio::Range->new(-start=>10,
                            -end=>20,
						    -strand=>1);

isa_ok($range,'Bio::Range', 'BioRange object');
is($range->strand, 1);

my $range2 = Bio::Range->new(-start=>15,
                             -end=>25,
						     -strand=>1);

isa_ok($range2,'Bio::Range', 'BioRange object');
is($range2->strand, 1);

my $r = Bio::Range->new();
is ( $r->strand(0), 0 ) ;
is ( $r->start(27), 27 );
is ( $r->end(28), 28 ) ;

ok(! defined $r->intersection($range2));

$r = $range->union($range2);
is($r->start, 10);
is($r->end, 25);

$r = $range->intersection($range2);
is ( $r->start, 15  ) ;
is ( $r->end, 20    );
is ( $r->strand, 1  );

# intersection and union can also take lists
my $range3 = Bio::Range->new(-start=>18,-end=>30);
isa_ok($range3,'Bio::Range', 'BioRange object');

$r = $range->intersection([$range2, $range3]);
ok( ( $r->start == 18 ) && ( $r->end == 20 ));
$r = Bio::Range->intersection([$range, $range2, $range3]);
ok($r->start == 18 && $r->end == 20);
$r = $range->union($range2, $range3);
ok( ( $r->start == 10 ) && ( $r->end == 30 ) );
$r = Bio::Range->union($range, $range2, $range3);
ok( ( $r->start == 10 ) && ( $r->end == 30 ) );
$range3->start(21);
ok (! $range->intersection([$range2, $range3]));

ok (! $range->contains($range2));
ok (! $range2->contains($range));
ok ($range->overlaps($range2));
ok ($range2->overlaps($range));

# testing strand
$range3 = Bio::Range->new(-start => 15,
                           -end => 25,
						   -strand => 1);

my $range4 = Bio::Range->new(-start => 15,
						    -end => 25,
							-strand => -1);

isa_ok($range4,'Bio::Range', 'BioRange object');

my $range5 = Bio::Range->new(-start => 15,
                             -end => 25,
						     -strand => 0);

isa_ok($range5,'Bio::Range', 'BioRange object');

my $range6 = Bio::Range->new(-start => 20,
							 -end => 30,
							 -strand => -1);

isa_ok($range6,'Bio::Range', 'BioRange object');

ok $range3->_ignore($range4), ' 1 & -1' ;     
ok $range3->_weak($range3),' 1 & 1 true' ;       
ok $range3->_weak($range5), ' 1 & 0 true' ;       
ok (! $range3->_weak($range4), ' 1 & -1 false' );   
ok $range3->_strong($range3), ' 1 & 1 true' ;     
ok (! $range3->_strong($range5), ' 1 & 0 false' ); 
ok (! $range3->_strong($range4), ' 1 & -1 false' ); 

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
is ( $r->strand, 0 );

# or if both strands are -1 then -1 should be returned
$r = $range6->intersection($range4);
is ( $r->strand, -1 );

# test implemention of offsetStranded:
$r = Bio::Range->new(-start => 30, -end => 40, -strand => -1);
isa_ok($r, 'Bio::Range', 'Bio::Range object') ;
is ($r->offsetStranded(-5,10)->toString, '(20, 45) strand=-1');
is ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=-1');
$r->strand(1);
is ($r->offsetStranded(-5,10)->toString, '(25, 50) strand=1');
is ($r->offsetStranded(+5,-10)->toString, '(30, 40) strand=1');
