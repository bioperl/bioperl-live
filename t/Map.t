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
    plan tests => 47;
}

END {
}

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
				      -position => [$map, 200] # THIS DOES NOT GET SET!
				      );
skip $marker2->map, $map;
#skip ($marker2->position->value, 200);

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

ok $marker3->known_maps, qw ( 1 ); 
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

