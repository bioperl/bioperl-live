# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG);
    $DEBUG = $ENV{'BIOPERLDEBUG'};
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 109;
}

END {
}

#
# Let's test first the map class : Bio::Map::CytoMap
#

use Bio::Map::CytoMap;
ok 1;

ok my $map = new Bio::Map::CytoMap(-name  => 'my');
ok $map->type, 'cyto'; 
ok $map->units, ''; 
ok $map->length, undef;
ok $map->name, 'my';
ok $map->species('human'), 'human';
ok $map->species, 'human';
ok $map->unique_id, '1';
#
#
# Secondly, we make sure the location calculations in
#           Bio::Map::CytoPosition make sense
#

use Bio::Map::CytoPosition;
use Bio::Range;

ok(1);

my($a, $b, $r);
my $string = 'b';
ok Bio::Map::CytoPosition::_pad($string, 5, 'z'), 'bzzzz';

ok $a = Bio::Map::CytoPosition->new();
ok $a->isa('Bio::Map::CytoPosition');
ok $a->cytorange, undef;


$a->verbose(2);
eval {
    ok $a->value('C'), 'C'; 
    ok $a->cytorange, undef ;
};
ok $@;
$a->verbose(0);

ok $a->value('X'), 'X'; 
$r = $a->cytorange;
ok $r->isa('Bio::Range');
ok $r->start, 100000000;
ok $r->end, 100200000;

$a->value('1p');
ok $a->cytorange->start, 1000000;
ok $a->cytorange->end, 1100000;

$a->value('2qter');
ok $a->cytorange->start, 2200000;
ok $a->cytorange->end, 2200000;

$a->value('2qcen');
ok $a->cytorange->start, 2100000;
ok $a->cytorange->end, 2100000;

eval {
    $a->value('2qcen2');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('2q22');
ok $a->cytorange->start, 2122000;
ok $a->cytorange->end, 2122999;

$a->value('2p22');
ok $a->cytorange->start, 2077001;
ok $a->cytorange->end, 2078000;

$a->value('2p21');
ok $a->cytorange->start, 2078001;
ok $a->cytorange->end, 2079000;

$a->value('10p22.1-cen');
ok $a->cytorange->start, 10022199;
ok $a->cytorange->end, 10100000;

eval {
    $a->value('10q22.1-cen');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('10q22.1-ter');
ok $a->cytorange->start, 10122100;
ok $a->cytorange->end, 10200000;


eval {
    $a->value('10q22.1-p');
    $a->cytorange->start;
};
ok 1 if $@;

$a->value('10qcen-qter');
ok $a->cytorange->start, 10100000;
ok $a->cytorange->end, 10200000;

$a->value('10pcen-qter');
ok $a->cytorange->start, 10100000;
ok $a->cytorange->end, 10200000;

$a->value('10q22.1-q23');
ok $a->cytorange->start, 10122100;
ok $a->cytorange->end, 10123999;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('10p22.1-p23');
ok $a->cytorange->start, 10076001;
ok $a->cytorange->end,  10077900;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('10cen-p23');
ok $a->cytorange->start, 10076001;
ok $a->cytorange->end, 10100000;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('10q22.1-p23');
ok $a->cytorange->start, 10076001;
ok $a->cytorange->end, 10122199;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('10p22.1-q23');
ok $a->cytorange->start, 10077801;
ok $a->cytorange->end, 10123999;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('10q22.1-p22');
ok $a->cytorange->start, 10077001 ;
ok $a->cytorange->end, 10122199 ;

$b = Bio::Map::CytoPosition->new();
$b->value('10p22-p22.1');
ok $b->cytorange->start, 10077801 ;
ok $b->cytorange->end, 10078000;
ok $a->cytorange->overlaps($b->cytorange);


$a->value('10p22.1-q23');
ok $a->cytorange->start, 10077801;
ok $a->cytorange->end, 10123999;
ok ($a->cytorange->start < $a->cytorange->end );


$a->value('17p13-pter');
ok $a->cytorange->start, 17000000;
ok $a->cytorange->end, 17087000;
ok ($a->cytorange->start < $a->cytorange->end );

$a->value('17cen-pter');
ok $a->cytorange->start, 17000000;
ok $a->cytorange->end, 17100000;
ok ($a->cytorange->start < $a->cytorange->end );


#-----------------------------------------
#my $s;

sub test {
    my ($s) = @_;
    my $a = Bio::Map::CytoPosition->new();
    $a->value($s);
    $r = $a->cytorange;
    ok $a->range2value($r), $s;
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


use Bio::Map::CytoMarker;
ok 1;

ok my $marker1 = new Bio::Map::CytoMarker();
ok $marker1->name('gene1'), 'Unnamed marker' ;
ok $marker1->name(), 'gene1';
ok $marker1->position($map, '10p33.13-q15');

ok my $marker2 = new Bio::Map::CytoMarker(-name => 'gene2' );
ok $marker2->position($map, '10p10-15');
ok $marker1->get_chr, 10;

ok my $marker3 = new Bio::Map::CytoMarker(-name => '3' );
ok $marker3->position($map, '10p1');

ok my $marker4 = new Bio::Map::CytoMarker(-name => '4' );
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

