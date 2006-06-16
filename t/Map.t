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
    plan tests => 74;
}

#END {
#}

#
# Original set of tests for 1.5.1
#

#
# Let's test first the map class : Bio::Map::SimpleMap
#
use Data::Dumper;
use Bio::Map::SimpleMap;
ok 1;

ok my $map = new Bio::Map::SimpleMap(-name  => 'my');
ok $map->type('cyto'); 
ok $map->type, 'cyto'; 
ok $map->units, ''; 
ok $map->length, 0, "Length is ". $map->length;
ok $map->name, 'my';
ok $map->species('human'), 'human';
ok $map->species, 'human';
ok $map->unique_id, '1';


#
# Secondly, create Markers for the Map
#

use Bio::Map::Marker;
ok 1;

#
# Four ways of adding a position:
#
# 1. position gets a map object and a value as arguments

ok my $marker1 = new Bio::Map::Marker();
ok $marker1->name('gene1'), 'Unnamed marker' ;
ok $marker1->name(), 'gene1';

ok $marker1->position($map, 100);
ok $marker1->position->value, 100;
ok $marker1->map, $map;


#
# Now that we have a Marker, let's 
#    make sure the basic Position class works
#

use Bio::Map::Position;
ok 1;
ok my $pos = new Bio::Map::Position();
ok $pos->map($map);
ok $pos->map(), $map;
ok $pos->marker($marker1);
ok $pos->marker(), $marker1;

ok $pos->value('999');
ok $pos->value(), '999';
ok $pos->numeric, 999;

# ... and then continue testing the Marker

# 2. position is set in the constructor
ok my $marker2 = new Bio::Map::Marker(-name => 'gene2',
				      -position => [$map, 200]
				      );
ok ( $marker2->map, $map);
ok ($marker2->position->value, 200);

# 3. marker is first added into map, 
#    then marker knows the the position belongs to
ok my $marker3 = new Bio::Map::Marker(-name => 'gene3');
ok $map->add_element($marker3);
ok $marker3->map, $map;
ok $marker3->position(300);
ok $marker3->position->value, 300;
ok my $position = $marker3->position($map);

# 4. A Position is explicitely created
ok my $cpos = new Bio::Map::Position(-map => $map,
				     -value => 500);
ok $marker3->position($cpos);
my $cpos2 = $marker3->position($cpos);
ok $cpos2 eq $cpos;
ok $marker3->position->value, 500;


#
# Next, what do markers know about Maps?
#

ok (scalar ($marker3->known_maps), 1); 
ok $marker3->in_map(1);
ok ! $marker3->in_map(2);



#
# Lastly, let's test the comparison methods
#

ok $marker1->equals($marker1);
ok ! $marker1->equals($marker3);
ok $marker1->less_than($marker3);
ok ! $marker1->greater_than($marker3);
ok ! $marker3->less_than($marker1);
ok $marker3->greater_than($marker1);

#
# 48 test so far...
# New set of tests for additions relating to complete implementation of everything suggested by the docs
#

#
# Let's try and do everything in the synopsis of Bio::Map::Marker
#
my $map1 = new Bio::Map::SimpleMap(-name => 'genethon', -type => 'Genetic');
my $pos1 =  new Bio::Map::Position(-map => $map1, -value => 100 );
my $pos2 = new Bio::Map::Position(-map => $map1, -value => 200 );
ok my $o_usat = new Bio::Map::Marker (-name=>'Chad Super Marker 2', -positions => [ [$map1, 100], [$map1, 200] ] );
ok my $o_usat2 = new Bio::Map::Marker (-name=>'Chad Super Marker 3', -positions => [ $pos1, $pos2 ] );
my $marker1 = new Bio::Map::Marker (-name=>'hypervariable1', -map => $map1, -position => 100 );
my $marker2 = new Bio::Map::Marker (-name=>'hypervariable2');
$map1->add_element($marker2);
ok $marker2->position(150);
ok $marker2->add_position(200);
my $map2 = new Bio::Map::SimpleMap(-name => 'xyz', -type => 'Genetic');
ok $marker2->add_position($map2,200);
ok $marker1->position->value, 100;
my ($pos1a, $pos2a) = $marker2->each_position($map1);
ok $pos1a->value == 150 && $pos2a->value == 200;

#
# 57 tests so far...
# Now test the methods in Bio::Map::Marker that haven't been tested yet
#

$marker2->positions([[$map1, 300], [$map2, 400]]);
my @p = map($_->numeric, $marker2->each_position);
ok $p[0] == 150 && $p[1] == 200 && $p[2] == 200 && $p[3] == 300 && $p[4] == 400;
$marker2->purge_positions($map2);
@p = map($_->numeric, $marker2->each_position);
ok $p[0] == 150 && $p[1] == 200 && $p[2] == 300;
@p = map($_->numeric, $marker1->less_than($marker2));
ok @p == 1 && $p[0] == 100;
@p = map($_->numeric, $marker2->greater_than($marker1));
ok @p == 3 && $p[0] == 150 && $p[1] == 200 && $p[2] == 300;

#
# 59 tests so far...
# Now test the new range-related things in position and marker
#

ok my $r_pos1 = new Bio::Map::Position(-map => $map1, -marker => $marker3, -start => 100, -length => 11);
ok $pos1->numeric, $r_pos1->numeric;
ok $pos1->start, $r_pos1->start;
ok $pos1->end, $r_pos1->start;
ok my $r_pos2 = new Bio::Map::Position(-map => $map1, -marker => $marker3, -start => 115, -end => 125);
my $r_pos3 = Bio::Map::Position->new(-map => $map, -start => 101, -end => 109);
my $r_pos4 = Bio::Map::Position->new(-map => $map, -start => 120, -end => 130);
my $r_pos5 = Bio::Map::Position->new(-map => $map, -start => 140, -end => 150);
my $marker3 = new Bio::Map::Marker (-name=>'qwe', -positions => [$r_pos1, $r_pos2]);
my $marker4 = new Bio::Map::Marker (-name=>'iop', -positions => [$r_pos3, $r_pos4, $r_pos5]);
ok @p = map($_->toString, $marker3->overlaps($marker4));
ok @p == 2 && $p[0] eq '100..110' && $p[1] eq '115..125';
ok @p = map($_->toString, $marker3->contains($marker4));
ok @p == 1 && $p[0] eq '100..110';
ok my $int = $r_pos2->intersection($r_pos4);
ok $int->toString eq '120..125';
my $r_pos6 = Bio::Map::Position->new(-map => $map, -start => 145, -end => 155);
my $r_pos7 = Bio::Map::Position->new(-map => $map, -start => 148, -end => 160);
ok $int = Bio::Map::Position->intersection([$r_pos5, $r_pos6, $r_pos7]);
ok $int->toString eq '148..150';
ok my $union = Bio::Map::Position->union($r_pos5, $r_pos6, $r_pos7);
ok $union->toString eq '140..160';
