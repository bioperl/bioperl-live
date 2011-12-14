# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(    #-tests => 10,
	-requires_module => 'LWP::UserAgent'
    );

    use_ok('Bio::DB::SeqVersion');
}

my $DEBUG = $ENV{BIOPERLDEBUG} || 0;

ok my $query = Bio::DB::SeqVersion->new( -type => 'gi' );

SKIP: {
    test_skip( -tests => 8, -requires_networking => 1 );

    throws_ok { $query->get_history('DODGY_ID_WHICH_SHOULD_FAIL') }
    qr/ID likely does not exist/i, 'throw on bad ID';

    #my $latest_gi = $query->get_recent(2);
    #is($latest_gi, 2, 'get_recent');
    #
    #my @all_gis = $query->get_all(2);
    #cmp_ok(@all_gis, '>=', 8, 'get_all');

    #my $latest_gi = $query->get_recent('A00002');
    #is($latest_gi, 2, 'get_recent, string');
    #
    #$latest_gi = $query->get_recent(27478738);
    #is($latest_gi, 42659163, 'get_recent, integer');
    #
    ## check that default type is "gi"
    #ok $query = Bio::DB::SeqVersion->new();
    #ok my $ref = $query->get_history(3245);
    #is($ref->[0]->[0], 578167, 'get_history');
}

done_testing();
