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
    plan tests => 7;
}

END {
}
#require 'dumpvar.pl';

use Bio::Map::LinkagePosition;
ok(1);

my $position = new Bio::Map::LinkagePosition(-positions => 2,
	-distance => 22.3);
ok(1);

ok(($position->each_position)[0], 2);
ok($position->distance, 22.3);

my $position2 = new Bio::Map::LinkagePosition(-positions => [qw(3 4 5)]);
ok(1);
# print("position2 looks like this:\n");
# dumpValue($position2);
ok($position2->distance() == 0);
ok($position2->each_position, 3);
