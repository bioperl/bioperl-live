# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {     
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}

	use Test::More;
	plan tests => 52; 
	use_ok('Bio::MapIO');
	use_ok('Bio::Root::IO');
}


my $verbose = $ENV{'BIOPERLDEBUG'};

ok my $mapio = new Bio::MapIO(-verbose => $verbose,
										-format => 'mapmaker',
										-file   => Bio::Root::IO->catfile
										('t','data','mapmaker.out'));

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

ok $mapio = new Bio::MapIO(-format => 'mapmaker',
									-file   => Bio::Root::IO->catfile
									('t','data','mapmaker.txt'));

$map = $mapio->next_map;
is $map->length,382.5;

$count = 0;
foreach my $marker ( $map->each_element ) {
	$count++;
	is($marker->position->order,$count);
}
is $count,13;

ok $mapio = Bio::MapIO->new(-format => 'fpc',
									-file   => Bio::Root::IO->catfile
									('t','data','ctgdemo.fpc'));

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
