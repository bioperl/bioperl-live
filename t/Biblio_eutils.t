# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $NUMTESTS = 5;
    plan tests => $NUMTESTS;
    eval { require 'IO/String.pm' };
    if( $@ ) {
	if( $DEBUG ) {
	    print STDERR "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n" if($DEBUG);
	}
	for( $Test::ntest..$NUMTESTS ) {
	    skip("IO::String not installed. Skipping tests",1);
	}
       $error = 1; 
    }
}

if( $error ==  1 ) {
    exit(0);
}
END{ 
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('unable to run all of the Biblio_biofetch tests',1);
    }
}
use Bio::Biblio;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $db;

my $verbose =  $DEBUG || 0;

eval { 
    ok ($db = new Bio::Biblio (-access => 'eutils',
			       -verbose=>$verbose));
    ok(defined($db->find('"Day A"[AU] AND ("Database Management Systems"[MH] OR "Databases, Genetic"[MH] OR "Software"[MH] OR "Software Design"[MH])')));
};

if ($@) {
    if( $DEBUG  ) { 
	print STDERR "Warning: Couldn't connect to Eutils server!\n" . $@;
    }
    foreach ( $Test::ntest..$NUMTESTS) { 
	skip('No network access - could not connect to PubMed Eutils',1);
    }
    exit(0);
}

while(my $xml = $db->get_next){
  ok(1);
}
