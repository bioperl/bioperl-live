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
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;
    
    plan tests => ($NUMTESTS = 5);
    eval { require IO::String; };
    if( $@ ) {
	warn( "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n") if $DEBUG;
    	$msg .= 'IO::String not installed. ';
	$error = 1;
    } 
    eval { require LWP::Simple; };
    if( $@ ) {
	warn( "LWP::Simple not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n") if $DEBUG;
	$msg .= 'LWP::Simple not installed. ';
	$error = 1; 
    }
}

exit(0) if $error;

END { 
    foreach ( $Test::ntest..$NUMTESTS) {
	skip($msg,1);
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
    warn "Warning: Couldn't connect to Eutils server!\n$@\n" if $DEBUG;
    $msg = 'No network access - could not connect to PubMed Eutils';
    exit(0);
}

while(my $xml = $db->get_next) {
    ok(1);
}
