# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Biblio.t'

use strict;
use vars qw($NUMTESTS);

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

    $NUMTESTS = 2;
    plan tests => $NUMTESTS;
    unless( eval "require SOAP::Lite; 1;" ) {
	print STDERR "SOAP::Lite not installed. This means that Bio::Biblio module is not usable. Skipping tests.\n";
	for( 1..$NUMTESTS ) {
	    skip (1,"SOAP::Lite not installed. This means that Bio::Biblio module is not usable. Skipping tests.\n");
	}
	$error = 1;
    }
}

if( $error ==  1 ) {
    exit(0);
}

require Bio::Biblio;

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

## !!!TBD: WORK IN PROGRESS!!!

# check 'use ...'
print "test 'use Bio::Biblio ...'\t";
eval 'use Bio::Biblio 99.99'; # hm, definitely should fail
ok ($@ =~ /99\.99 required/);

# check 'new...'
my $biblio;
print "test 'new Bio::Biblio ...'\t";
ok defined ($biblio = new Bio::Biblio (-location => 'http://localhost:4567'));


