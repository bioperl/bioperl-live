# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $error);


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

    $NUMTESTS = 1;
    plan tests => $NUMTESTS;
    eval { require IO::String; 
	   require Bio::Tools::Phylo::PAML;}; 
    if( $@ ) { print STDERR "no IO string installed\n"; 
	$error = 1;
	}
}

END { 
    foreach ( $Test::ntest .. $NUMTESTS ) {
	skip("unable to run all of the PAML tests",1);
    }
}


exit(0) if( $error );

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $inpaml = new Bio::Tools::Phylo::PAML(-file => 't/data/codeml.mlc');
ok($inpaml);
