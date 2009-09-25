# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 110);
    
    use_ok('Bio::Map::CytoMap');
    use_ok('Bio::Map::CytoPosition');
    use_ok('Bio::Map::CytoMarker');
}

#
# Let's test first the map class : Bio::Map::CytoMap
#

ok my $map = Bio::Map::CytoMap->new(-name  => 'my');
is $map->type, 'cyto'; 
is $map->units, ''; 
is $map->length, 0;
is $map->name, 'my';
is $map->species('human'), 'human';
is $map->species, 'human';
is $map->unique_id, '1';

#
#
# Secondly, we make sure the location calculations in
#           Bio::Map::CytoPosition make sense
#

my($a, $b, $r);
my $string = 'b';
is Bio::Map::CytoPosition::_pad($string, 5, 'z'), 'bzzzz';

ok $a = Bio::Map::CytoPosition->new();
isa_ok $a, 'Bio::Map::CytoPosition';
is $a->cytorange, undef;


$a->verbose(2);
eval {
    is $a->value('C'), 'C'; 
    is $a->cytorange, undef ;
};
ok $@;
$a->verbose(0);

is $a->value('X'), 'X';
$r = $a->cytorange;
isa_ok $r, 'Bio::Range';
is $r->start, 100000000;
is $r->end, 100200000;

$a->value('1p');
is $a->cytorange->start, 1000000;
is $a->cytorange->end, 1100000;

$a->value('2qter');
is $a->cytorange->start, 2200000;
is $a->cytorange->end, 2200000;

$a->value('2qcen');
is $a->cytorange->start, 2100000;
is $a->cytorange->end, 2100000;

eval {
    $a->value('2qcen2');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('2q22');
is $a->cytorange->start, 2122000;
is $a->cytorange->end, 2122999;

$a->value('2p22');
is $a->cytorange->start, 2077001;
is $a->cytorange->end, 2078000;

$a->value('2p21');
is $a->cytorange->start, 2078001;
is $a->cytorange->end, 2079000;

$a->value('10p22.1-cen');
is $a->cytorange->start, 10022199;
is $a->cytorange->end, 10100000;

eval {
    $a->value('10q22.1-cen');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('10q22.1-ter');
is $a->cytorange->start, 10122100;
is $a->cytorange->end, 10200000;


eval {
    $a->value('10q22.1-p');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('10qcen-qter');
is $a->cytorange->start, 10100000;
is $a->cytorange->end, 10200000;

$a->value('10pcen-qter');
is $a->cytorange->start, 10100000;
is $a->cytorange->end, 10200000;

$a->value('10q22.1-q23');
is $a->cytorange->start, 10122100;
is $a->cytorange->end, 10123999;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('10p22.1-p23');
is $a->cytorange->start, 10076001;
is $a->cytorange->end,  10077900;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('10cen-p23');
is $a->cytorange->start, 10076001;
is $a->cytorange->end, 10100000;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('10q22.1-p23');
is $a->cytorange->start, 10076001;
is $a->cytorange->end, 10122199;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('10p22.1-q23');
is $a->cytorange->start, 10077801;
is $a->cytorange->end, 10123999;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('10q22.1-p22');
is $a->cytorange->start, 10077001 ;
is $a->cytorange->end, 10122199 ;

$b = Bio::Map::CytoPosition->new();
$b->value('10p22-p22.1');
is $b->cytorange->start, 10077801 ;
is $b->cytorange->end, 10078000;
ok $a->cytorange->overlaps($b->cytorange);

$a->value('10p22.1-q23');
is $a->cytorange->start, 10077801;
is $a->cytorange->end, 10123999;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('17p13-pter');
is $a->cytorange->start, 17000000;
is $a->cytorange->end, 17087000;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );

$a->value('17cen-pter');
is $a->cytorange->start, 17000000;
is $a->cytorange->end, 17100000;
cmp_ok ($a->cytorange->start, '<', $a->cytorange->end );


#-----------------------------------------
#my $s;

sub test {
    my ($s) = @_;
    my $a = Bio::Map::CytoPosition->new();
    $a->value($s);
    $r = $a->cytorange;
    is $a->range2value($r), $s;
}

test '1';
test '2p';
test '3q';
test '4cen';

test '5pter';
test '6qter';
test '7p21';
test '8q11.1';

test '9q13.13-15';
test '10p13.13-q15';
test '11p13.13-qter';
test '12p13.13-qter';
test '13p13.13-14';
test '14p13.13-pter';
test '15cen-q2';
test '16cen-p2';
#test '17cen-pter'; eq 17p
#test '18cen-qter'; eq 18q 


# by now we should be convinced that band conversion to a range works
# so lets try to use it for comparing markers.

ok my $marker1 = Bio::Map::CytoMarker->new();
is $marker1->name('gene1'), 'gene1' ;
ok $marker1->position($map, '10p33.13-q15');

ok my $marker2 = Bio::Map::CytoMarker->new(-name => 'gene2' );
ok $marker2->position($map, '10p10-15');
is $marker1->get_chr, 10;

ok my $marker3 = Bio::Map::CytoMarker->new(-name => '3' );
ok $marker3->position($map, '10p1');

ok my $marker4 = Bio::Map::CytoMarker->new(-name => '4' );
ok $marker4->position($map, '10q2');

#
# Lastly, let's test the comparison methods
#

ok $marker1->equals($marker1);
ok ! $marker1->equals($marker2);

ok $marker3->less_than($marker4);
ok ! $marker3->greater_than($marker4);
ok ! $marker4->less_than($marker3);
ok $marker4->greater_than($marker3);

ok ! $marker4->overlaps($marker3);
ok $marker1->overlaps($marker3);

ok ! $marker4->contains($marker3);
ok $marker1->contains($marker3);

# 
# Test throw() in some private functions
#

eval { Bio::Map::CytoPosition::_pad('string', -1, 'x'); };
like($@, qr/positive integer/);
eval { Bio::Map::CytoPosition::_pad('string', +1, 'toolong'); };
like($@, qr/single character/);

