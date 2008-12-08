# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 51);
	
	use_ok('Bio::MapIO');
}

my $verbose = test_debug();

ok my $mapio = Bio::MapIO->new(-verbose => $verbose,
										-format => 'mapmaker',
										-file   => test_input_file('mapmaker.out'));

my $map = $mapio->next_map;

isa_ok($map, 'Bio::Map::MapI');

is $map->units, 'cM';
is $map->type, 'Genetic';
is $map->name('test map'), 'test map'; # map name is unset for this data type

my $count;
foreach my $marker ( $map->each_element ) {
	$count++;
	is($marker->position->order,$count);
}
is $count,18;

ok $mapio = Bio::MapIO->new(-format => 'mapmaker',
									-file   => test_input_file('mapmaker.txt'));

$map = $mapio->next_map;
is $map->length,382.5;

$count = 0;
foreach my $marker ( $map->each_element ) {
	$count++;
	is($marker->position->order,$count);
}
is $count,13;

ok $mapio = Bio::MapIO->new(-format => 'fpc',
									-file   => test_input_file('ctgdemo.fpc'));

$map = $mapio->next_map;

is($map->length, 0);
is($map->version, 7.2);
is($map->modification_user, 'cari');
is($map->group_type, 'Chromosome');
is($map->group_abbr, 'Chr');
is($map->core_exists, 0);

$count = 0;
foreach my $marker ($map->each_markerid) {
	$count++;
}

is($count,150);

# add tests for get_markerobj

$count = 0;
foreach my $clone ($map->each_cloneid) {
	$count++;
}

is($count,618);

# add tests for get_cloneobj

$count = 0;
foreach my $contig ($map->each_contigid) {
	$count++;
}

is($count,2);

# add tests for get_contigobj

# need tests for
# matching_bands
# coincidence_score
# print_contiglist
# print_markerlist
# print_gffstyle
