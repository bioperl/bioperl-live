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
    eval { require Test::More; };
    $error = $DEBUG ? '' : 'Must set BIOPERLDEBUG=1 for network tests';
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    eval { require IO::String; };
	if( $@ ) {
		$error = "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.";
	} 
	eval { require LWP::Simple; };
	if( $@ ) {
		$error = "LWP::Simple not installed. This means the Bio::DB::* modules are not usable. Skipping tests.";
	}
	
	if ($error) {
		plan skip_all => $error;
	} else {
		plan tests => ($NUMTESTS = 6);
	}
	use_ok('Bio::Biblio');
}

## End of black magic.

my $db;

my $verbose =  $DEBUG || 0;
SKIP: {
	ok ($db = Bio::Biblio->new(-access => 'eutils',
					   -verbose=>$verbose));
	eval { 
		ok(defined($db->find('"Day A"[AU] AND ("Database Management Systems"[MH] OR "Databases, Genetic"[MH] OR "Software"[MH] OR "Software Design"[MH])')));
	};
	
	if ($@) {
		skip("Warning: Couldn't connect to Eutils server!\n$@\n",1);
	}
	
	# these aren't exactly the most stringent of tests...
	
	while(my $xml = $db->get_next) {
		ok(1);
	}
}