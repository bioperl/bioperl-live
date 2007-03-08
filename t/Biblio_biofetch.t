# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG $error $msg);

BEGIN { 
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test::More; };
	$error = 0;
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
 
	plan tests => ($NUMTESTS = 13);

	use_ok('Bio::Biblio');
	use_ok('Bio::Biblio::IO');
}

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my ($db,$ref,$refio);
# get a single ref

my $verbose =  $DEBUG || 0;

$ref = $refio = undef;

# check BioFetch access method
SKIP: {
	eval { require IO::String; };
	if( $@ ) {
		skip( "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.",11);
	}
	eval { require LWP::Simple; };
	if( $@ ) {
		skip( "LWP::Simple not installed. This means the Bio::DB::* modules are not usable. Skipping tests.",11);
	}
	eval { require HTTP::Request::Common; };
	if( $@ ) {
		skip( "HTTP::Request::Common not installed. This means the Bio::DB::* modules are not usable. Skipping tests.",11);
	}
	ok ($db = new Bio::Biblio (-access => 'biofetch',
										# -verbose => $verbose,
									  ));
	eval { 
		$ref = $db->get_by_id('10592273');
	};

	if ($@) {
		skip( "Warning: Couldn't connect to BioFetch server with Bio::DB::Biblio::biofetch!$@",10); 
	}
	ok(defined($ref)); 
	is $ref->identifier, '10592273';
	$ref = $refio = undef;

	ok defined($db = new Bio::Biblio(-access => 'biofetch',
												# -verbose => $verbose,
											   )); 

	my $ids = ['10592273', '9613206'];
	eval {
		$refio = $db->get_all($ids);
	};
	
	if ($@) {
		skip("Batch access test failed.Error: $@",7);
	}
	
	ok(defined($refio));
	is($refio->next_bibref->identifier, '9613206');
	is($refio->next_bibref->identifier, '10592273');

	ok defined($db = new Bio::Biblio(-access => 'biofetch',
												# -verbose => $verbose,
											  )); 
	eval {
		$refio = $db->get_Stream_by_id(['10592273', '9613206']);
	};

	if ($@) {    
		skip("Batch access test failed.Error: $@",3);
	}
	
	ok(defined($refio));	
	is($refio->next_bibref->identifier, '9613206');
	is($refio->next_bibref->identifier, '10592273');
}
