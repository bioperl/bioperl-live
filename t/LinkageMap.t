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
    plan tests => 16;
}

END {
}
require 'dumpvar.pl';

use Bio::Map::LinkagePosition;
use Bio::Map::Microsatellite;
use Bio::Map::LinkageMap;

ok 1 ;
my $verbose = 0;
ok my $map = new Bio::Map::LinkageMap('-verbose' => $verbose,
				   '-name'    => 'Leviathon',
				   '-type'    => 'Genetic',
				   '-units'   => 'cM',
				   '-species' => "Brassica");
ok ref($map), 'Bio::Map::LinkageMap';
ok $map->name, 'Leviathon';
ok $map->type, 'Genetic';
ok $map->units, 'cM';
ok $map->species, 'Brassica';
ok $map->unique_id, '1';

ok my $position = new Bio::Map::LinkagePosition('-order' => 2,
						'-map' =>  $map, 
						'-value' => 22.3
						);

ok $position->order, 2;
ok $position->map, $map,
ok $position->value, 22.3;

ok my $o_usat = new Bio::Map::Microsatellite('-name'     => "Chad marker",
					     '-position' => $position);

ok $o_usat->name, 'Chad marker';
ok $o_usat->position, $position ;
ok $map->add_element($o_usat);

#use Data::Dumper; print Dumper($map);
#----------------------------
#ok my $position2 = new Bio::Map::LinkagePosition(-order => qw(3 4 5),
#						 );
# print("position2 looks like this:\n");
# dumpValue($position2);
#ok(($position2->each_position_value('fakemap'))[0] == 0);
#ok $position2->order, 3;

#-------------
#ok($position->order, 2);
#ok(($position->each_position_value($map))[0], 22.3);
	# what should be printed if this was ok?
	# ok(1);

#ok my $o_usat = new Bio::Map::Microsatellite('-name'     => "Chad marker",
#					      '-position' => $position);
#
#ok $o_usat->name, 'Chad marker';
#ok $o_usat->position, $position ;
#ok $map->add_element($o_usat);
# what should be printed if this is ok?
#dumpValue($map);

# add more tests
# see also t/microsatellite.t and t/linkageposition.t
