# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($DEBUG $TESTCOUNT);
$DEBUG = $ENV{'BIOPERLDEBUG'};

BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	$TESTCOUNT = 6;
	plan tests => $TESTCOUNT;
}

use Bio::DB::SeqVersion;

ok 1;

if ($DEBUG) {
	my $query = Bio::DB::SeqVersion->new(-type => 'gi');

	my $latest_gi = $query->get_recent(2);
	ok($latest_gi,2);

	my @all_gis = $query->get_all(2);
	ok(scalar @all_gis,8);

	$latest_gi = $query->get_recent('A00002');
	ok($latest_gi,2);

	$latest_gi = $query->get_recent(27478738);
	ok($latest_gi,42659163);

	# check that default type is "gi"
	$query = Bio::DB::SeqVersion->new();
	my $ref = $query->get_history(3245);
	ok($ref->[0]->[0],578167);
} else {
	for ( $Test::ntest..$TESTCOUNT) {
		skip("Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test", 1);
        }
}
