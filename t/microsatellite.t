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

print("Checking if the Bio::Map::Microsatellite module could be used.\n");
use Bio::Map::Microsatellite;
ok(1);

print("Creating a microsatellite...\n");

my $o_usat = new Bio::Map::Microsatellite(-name=>'Chad Super Marker 2',
        -sequence => 'gctgactgatcatatatatatatatatatatatatatatatcgcgatcgtgatttt',
        -motif => 'at',
        -repeats => 15,
        -repeat_start_position => 12
        );

print("Getting the leading flanking sequence...\n");
ok($o_usat->get_leading_flank() eq "gctgactgatc");
print("Getting the trailing flanking sequence...\n");
ok($o_usat->get_trailing_flank() eq "cgcgatcgtgatttt");


dumpValue($o_usat);

