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
    
    plan tests => ($NUMTESTS = 11);
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
    eval { require HTTP::Request::Common; };
    if( $@ ) {
	warn( "HTTP::Request::Common not installed. This means the Bio::DB::* modules are not usable. Skipping tests.\n") if $DEBUG;
	$msg .= 'HTTP::Request::Common not installed. ';
	$error = 1; 
    }
}

END{ 
    foreach ( $Test::ntest..$NUMTESTS) {
	skip($msg,1);
    }
}

exit if $error;

use Bio::Biblio;
use Bio::Biblio::IO;

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


eval { 
    ok ($db = new Bio::Biblio (-access => 'biofetch',
			       -verbose=>$verbose));
    ok(defined($ref = $db->get_by_id('10592273')));
    ok $ref->identifier, '10592273';
};

if ($@) {

    warn( "Warning: Couldn't connect to BioFetch server with Bio::DB::Medline!\n$@\n") if $DEBUG;
    $msg = "Couldn't connect to BioFetch server with Bio::DB::Medline";
    exit(0);
}

$ref = $refio = undef;

eval {
    ok defined($db = new Bio::Biblio(-access => 'biofetch',
				     -verbose=>$verbose,
				     -retrievaltype => 'tempfile'
				    )); 


    my $ids = ['10592273', '9613206'];
    ok(defined($refio = $db->get_all($ids)));
    ok($refio->next_bibref->identifier, '10592273');
    ok($refio->next_bibref->identifier, '9613206');
};

if ($@) {    
    warn "Batch access test failed.Error: $@\n" if $DEBUG;
    $msg = 'No network access';
    exit(0);
}

eval {
    ok defined($db = new Bio::Biblio(-access => 'biofetch',
				     -verbose=>$verbose
				     )); 

    ok(defined($refio = $db->get_Stream_by_id(['10592273', '9613206'])));
    ok($refio->next_bibref->identifier, '10592273');
    ok($refio->next_bibref->identifier, '9613206');
};

if ($@) {    
    warn "Batch access test failed.Error: $@\n" if $DEBUG;
    $msg = 'No network access';
    exit(0);
}


