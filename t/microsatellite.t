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
    plan tests => 6;
}

END {
}
require 'dumpvar.pl';
use Bio::Map::SimpleMap;
use Bio::Map::Position;
use Bio::Map::Microsatellite;
ok(1);

my $map = new Bio::Map::SimpleMap(-units => 'MB',
				  -type  => 'oo-121');
my $position = new Bio::Map::Position(-positions => [ [$map, 20] ] );

my $o_usat = new Bio::Map::Microsatellite
    (-name=>'Chad Super Marker 2',
     -sequence => 'gctgactgatcatatatatatatatatatatatatatatatcgcgatcgtgatttt',
     -motif => 'at',
     -repeats => 15,
     -repeat_start_position => 12,
     -position => $position,
     );

ok($o_usat->get_leading_flank(), "gctgactgatc");
ok($o_usat->get_trailing_flank(), "cgcgatcgtgatttt");
ok($o_usat->motif(), 'at');
ok($o_usat->repeats(), 15);
ok($o_usat->repeat_start_position, 12);


#dumpValue($o_usat);

