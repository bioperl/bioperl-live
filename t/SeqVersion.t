# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($DEBUG $NUMTESTS);

BEGIN {
  $NUMTESTS = 10;
  $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

  eval { require Test::More; };
  if ($@) {
    use lib 't/lib';
  }
  use Test::More;
  
  eval {
	require LWP::UserAgent;
  };
  if ($@) {
	plan skip_all => 'LWP::UserAgent not installed. This means that the module is not usable. Skipping tests';
  }
  else {
	plan tests => $NUMTESTS;
  }
  
  use_ok('Bio::DB::SeqVersion');
}

ok my $query = Bio::DB::SeqVersion->new(-type => 'gi');

SKIP: {
	skip "Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test", 8 unless $DEBUG;

        eval { $query->get_history('DODGY_ID_WHICH_SHOULD_FAIL') };
        like($@, qr/could not parse/i, 'throw on bad ID');

	my $latest_gi = $query->get_recent(2);
	is($latest_gi, 2, 'get_recent');

	my @all_gis = $query->get_all(2);
	is(scalar @all_gis, 8, 'get_all');

	$latest_gi = $query->get_recent('A00002');
	is($latest_gi, 2, 'get_recent, string');

	$latest_gi = $query->get_recent(27478738);
	is($latest_gi, 42659163, 'get_recent, integer');

	# check that default type is "gi"
	ok $query = Bio::DB::SeqVersion->new();
	ok my $ref = $query->get_history(3245);
	is($ref->[0]->[0], 578167, 'get_history');
} 

