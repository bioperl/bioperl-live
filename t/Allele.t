# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

## We start with some black magic to print on failure.

use Test;
use strict;

BEGIN { plan tests => 15 }

use Bio::Variation::Allele;

ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my($a,$trunc,$rev);

$a = Bio::Variation::Allele->new(-seq=>'ACTGACTGACTG',
			-display_id => 'new-id',
			-moltype => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
ok defined $a,
ok ref($a), 'Bio::Variation::Allele';

ok $a->accession_number(), 'X677667';
ok $a->seq(), 'ACTGACTGACTG';
ok $a->display_id(),'new-id' ;
ok $a->desc, 'Sample Bio::Seq object';
ok $a->moltype(), 'dna';

ok defined($trunc = $a->trunc(1,4));
ok $trunc->seq(), 'ACTG', "Expecting ACTG. Got ". $trunc->seq();

ok defined($rev = $a->revcom());
ok $rev->seq(), 'CAGTCAGTCAGT';

$a->is_reference(1);
ok $a->is_reference;

$a->repeat_unit('ACTG');
ok $a->repeat_unit, 'ACTG';

$a->repeat_count(3);
ok $a->repeat_count, 3;

