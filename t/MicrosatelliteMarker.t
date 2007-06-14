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
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 9;
	use_ok('Bio::Map::SimpleMap');
	use_ok('Bio::Map::Position');
	use_ok('Bio::Map::Microsatellite');
}

require_ok 'dumpvar.pl';

my $map = Bio::Map::SimpleMap->new(-units => 'MB',
				  -type  => 'oo-121');
my $position = Bio::Map::Position->new(-map => $map,
				      -value => 20
				      );

my $o_usat = Bio::Map::Microsatellite->new
    (-name=>'Chad Super Marker 2',
     -sequence => 'gctgactgatcatatatatatatatatatatatatatatatcgcgatcgtgatttt',
     -motif => 'at',
     -repeats => 15,
     -repeat_start_position => 12,
     -position => $position,
     );

is($o_usat->get_leading_flank(), "gctgactgatc");
is($o_usat->get_trailing_flank(), "cgcgatcgtgatttt");
is($o_usat->motif(), 'at');
is($o_usat->repeats(), 15);
is($o_usat->repeat_start_position, 12);


#dumpValue($o_usat);

