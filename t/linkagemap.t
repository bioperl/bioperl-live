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
require 'dumpvar.pl';

use Bio::Map::LinkagePosition;
use Bio::Map::Microsatellite;

print("Checking if the Bio::Map::LinkageMap module could be used.\n");
        # test 1
use Bio::Map::LinkageMap;
ok(1);

print("Creating a new LinkageMap...\n");
my $map = new Bio::Map::LinkageMap(-name => 'Leviathon',
				-type => 'Genetic',
				-units=> 'cM',
				-species => "Brassica");
ok ref($map) eq 'Bio::Map::LinkageMap';


print("Creating a Position with a scalar...\n");
my $position = new Bio::Map::LinkagePosition(-positions => 2,
	-distance => 22.3);
	# what should be printed if this was ok?
	# ok(1);

print("Creating a Microsatellite marker with that position...\n");
my $o_usat = new Bio::Map::Microsatellite(-name=> "Chad marker",
	-position => $position);
	# what should be printed if this was ok?
	# ok(1);

print("Adding that to the LinkageMap...\n");
$map->add_element($o_usat);
	# what should be printed if this is ok?

# add more tests
# see also t/microsatellite.t and t/linkageposition.t
