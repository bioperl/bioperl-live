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
    plan tests => 22;
}

END {
}

use Bio::Map::LinkageMap;
use Bio::Map::Microsatellite;
use Bio::Map::LinkagePosition;
use Bio::Root::IO;
require 'dumpvar.pl';

use Bio::Map::Mapmaker;
ok(1);
my $verbose = -1;
my $in_mm = Bio::Map::Mapmaker->new
    ('-file' => Bio::Root::IO->catfile('t', 'data','mapmaker.out'),
     '-verbose' => $verbose
     );
ok(1);

my $summary_info;
eval {
    $summary_info = $in_mm->summary_info(); # this SHOULD fail
};
#print @$;

my $map = new Bio::Map::LinkageMap
   ('-verbose' => $verbose,
    '-name'    => "Chad Map",
    '-units'   => 'cM');

my ($marker,$marker_number,$marker_name,$marker_distance,$o_marker);
my $count = 1;
while ($marker = $in_mm->next_marker()) {
    ok(($marker->position->each_position)[0],$count++);
    # what kind of marker was that again?
#    bless $marker,"Bio::Map::Microsatellite";
    $map->add_element($marker);	
}

#dumpValue($map);

$summary_info = $in_mm->summary_info();
ok($map->name, "Chad Map");
ok($map->length, 107.8); # this is odd why is map length different for summary info
ok($map->units, 'cM');
print("The summary info is $summary_info\n");
#print("The total distance of the map ".$map->name()." is ".$map->length()." ".$map->units()."\n");

