# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

# Note this uses Test::More; this should catch the few perl versions w/o
# this test suite

use strict;
use lib '..','.','./blib/lib';
use vars qw($NUMTESTS $DEBUG $error);

BEGIN { 
	$NUMTESTS = 15;
	$error = 0;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0; 

    eval {require Test::More;};

	if ($@) {
		print STDERR 'Test::More is required for EUtilities.t, skipping...\n';
		use lib 't';
		use Test;
		$error = 1;
		foreach ( 1..$NUMTESTS) {
		   skip('Unable to run all of the DB tests',1);
		}
		exit(0);
	}
	
	if (!$DEBUG) {
		plan skip_all => 'Set BIOPERLDEBUG = 1 to run EUtilities.t tests';
	} else {
		plan tests => $NUMTESTS;
	}
}

END {
	if ($error) {
		foreach ( 1..$NUMTESTS) {
		   skip('Unable to run all of the DB tests',1);
		}
	}
}

exit(0) if $error;

use_ok('Bio::DB::EUtilities');
require_ok('LWP::UserAgent');
require_ok('XML::Simple');

# protein acc
my @acc = qw(MUSIGHBA1 P18584 CH402638);
# protein GI
my @ids = (70733119,68347418,26991676,82736201);

# EFetch

# ESearch

# ESearch (cookies)

# EInfo

# EPost 

# ELink (normal)

# ELink (multi_id)



# test 

#
# 
#

SKIP: {
	my $epost = Bio::DB::EUtilities->new(
                                    -eutil      => 'epost',
                                    -db         => 'protein',
                                    -id         => \@ids,
                                      );
		  
	isa_ok($epost, 'Bio::DB::EUtilities::epost');
	my $response;
	eval {$response = $epost->get_response; };
	skip "EPost HTTP error, not BioPerl's fault", 11 if $@;
	isa_ok($response, 'HTTP::Response');
	my $cookie = $epost->next_cookie;
	isa_ok($cookie, 'Bio::DB::EUtilities::Cookie');
	
	is($cookie->eutil, 'epost', 'cookie methods');
	is($cookie->database, 'protein', 'cookie methods');
	# these are not set using epost
	is($cookie->elink_dbfrom, undef, 'cookie methods');
	is($cookie->esearch_total, undef, 'cookie methods');
	is($cookie->esearch_query, undef, 'cookie methods');
	is($cookie->elink_queryids, undef, 'cookie methods');
	is($cookie->elink_linkname, undef, 'cookie methods');
	my ($webenv, $key) = @{ $cookie->cookie };
	like($webenv, qr{^\S{50}}, 'cookie WebEnv');
	like($key, qr{^\d+}, 'cookie query key');
}