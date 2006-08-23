# -*-Perl-*-
# Bioperl Test Harness Script for Modules
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
	plan tests => 51; 
}

if( $error == 1 ) {
	exit(0);
}

use Bio::MapIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'};
ok 1;

ok my $mapio = new Bio::MapIO(-verbose => $verbose,
										-format => 'mapmaker',
										-file   => Bio::Root::IO->catfile
										('t','data','mapmaker.out'));

my $map = $mapio->next_map;

ok(ref($map) && $map->isa('Bio::Map::MapI'));

ok $map->units, 'cM';
ok $map->type, 'Genetic';
ok $map->name('test map'), 'test map'; # map name is unset for this data type

my $count;
foreach my $marker ( $map->each_element ) {
	$count++;
	ok($marker->position->order,$count);
}
ok $count,18;

ok $mapio = new Bio::MapIO(-format => 'mapmaker',
									-file   => Bio::Root::IO->catfile
									('t','data','mapmaker.txt'));

$map = $mapio->next_map;
ok $map->length,382.5;

$count = 0;
foreach my $marker ( $map->each_element ) {
	$count++;
	ok($marker->position->order,$count);
}
ok $count,13;

ok $mapio = Bio::MapIO->new(-format => 'fpc',
									-file   => Bio::Root::IO->catfile
									('t','data','ctgdemo.fpc'));

$map = $mapio->next_map;

ok($map->length, 0);
ok($map->version, 7.2);
ok($map->modification_user, 'cari');
ok($map->group_type, 'Chromosome');
ok($map->group_abbr, 'Chr');
ok($map->core_exists, 0);

$count = 0;
foreach my $marker ($map->each_markerid) {
	$count++;
}

ok($count,150);

# add tests for get_markerobj

$count = 0;
foreach my $clone ($map->each_cloneid) {
	$count++;
}

ok($count,618);

# add tests for get_cloneobj

$count = 0;
foreach my $contig ($map->each_contigid) {
	$count++;
}

ok($count,2);

# add tests for get_contigobj

# need tests for
# matching_bands
# coincidence_score
# print_contiglist
# print_markerlist
# print_gffstyle
