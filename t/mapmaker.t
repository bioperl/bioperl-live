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
    plan tests => 12;
}

END {
}

use Bio::Map::LinkageMap;
use Bio::Map::Microsatellite;
use Bio::Map::LinkagePosition;

require 'dumpvar.pl';

print("Checking if the Bio::Map::Mapmaker module could be used.\n");
use Bio::Map::Mapmaker;
ok(1);

print("Checking to see if a mapmaker file can be parsed...\n");

my $in_mm = Bio::Map::Mapmaker->new('-file' => "t/data/mapmaker.out");
ok(1);

my $summary_info;
print("Checking to see if summary info can be pulled from the file before the markers...\n");
eval {
	$summary_info = $in_mm->summary_info();
};
print @$;

print("Creating a LinkageMap to place the markers on..\n");
my $map = new Bio::Map::LinkageMap(-name => "Chad's Superterriffic Linkage Map",
				-units => 'cM');

print("Checking to see if markers can be pulled from the file.\n");
my ($marker,$marker_number,$marker_name,$marker_distance,$o_marker);
while ($marker = $in_mm->next_marker()) {
		# what kind of marker was that again?
	bless $marker,"Bio::Map::Microsatellite";
	$map->add_element($marker);	
}

dumpValue($map);

$summary_info = $in_mm->summary_info();
print("The summary info is $summary_info\n");
print("The total distance of the map ".$map->name()." is ".$map->length()." ".$map->units()."\n");

