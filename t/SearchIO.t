# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 11; 
}

use Bio::SearchIO;

ok(1);

# test with RPSBLAST data first 
my $searchio = new Bio::SearchIO('-format' => 'blastxml',
				 '-file'   => 't/data/ecoli_domains.rps.xml');

while( my $report = $searchio->next_report )
{
    ok($report);
}


