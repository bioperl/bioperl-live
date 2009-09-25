# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
	use_ok('Bio::Map::LinkagePosition');
	use_ok('Bio::Map::Microsatellite');
	use_ok('Bio::Map::LinkageMap');
}

#require_ok('dumpvar.pl');

my $verbose = test_debug();
ok my $map = Bio::Map::LinkageMap->new('-verbose' => $verbose,
				   '-name'    => 'Leviathon',
				   '-type'    => 'Genetic',
				   '-units'   => 'cM',
				   '-species' => "Brassica");
is ref($map), 'Bio::Map::LinkageMap';
is $map->name, 'Leviathon';
is $map->type, 'Genetic';
is $map->units, 'cM';
is $map->species, 'Brassica';
is $map->unique_id, '1';

ok my $position = Bio::Map::LinkagePosition->new('-order' => 2,
						'-map' =>  $map, 
						'-value' => 22.3
						);

is $position->order, 2;
ok my $map2 = $position->map;
is $position->value, 22.3;

ok my $o_usat = Bio::Map::Microsatellite->new('-name'     => "Chad marker",
					     '-position' => $position);

is $o_usat->name, 'Chad marker';
is $o_usat->position, $position ;
ok $map->add_element($o_usat);

#use Data::Dumper; print Dumper($map);
#----------------------------
#ok my $position2 = Bio::Map::LinkagePosition->new(-order => qw(3 4 5),
#						 );
# print("position2 looks like this:\n");
# dumpValue($position2);
#is(($position2->each_position_value('fakemap'))[0], 0);
#is $position2->order, 3;

#-------------
#is($position->order, 2);
#is(($position->each_position_value($map))[0], 22.3);
	# what should be printed if this was ok?
	# ok(1);

#ok my $o_usat = Bio::Map::Microsatellite->new('-name'     => "Chad marker",
#					      '-position' => $position);
#
#is $o_usat->name, 'Chad marker';
#is $o_usat->position, $position ;
#ok $map->add_element($o_usat);
# what should be printed if this is ok?
#dumpValue($map);

# add more tests
# see also t/microsatellite.t and t/linkageposition.t
