# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) { 
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 14;
	use_ok('Bio::Variation::Allele');	
	}

my($a,$trunc,$rev);

$a = Bio::Variation::Allele->new(-seq=>'ACTGACTGACTG',
			-display_id => 'new-id',
			-alphabet => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
isa_ok($a, 'Bio::Variation::Allele');

is $a->accession_number(), 'X677667';
is $a->seq(), 'ACTGACTGACTG';
is $a->display_id(),'new-id' ;
is $a->desc, 'Sample Bio::Seq object';
is $a->alphabet(), 'dna';

ok defined($trunc = $a->trunc(1,4));
is $trunc->seq(), 'ACTG', $trunc->seq();

ok defined($rev = $a->revcom());
is $rev->seq(), 'CAGTCAGTCAGT';

$a->is_reference(1);
ok $a->is_reference;

$a->repeat_unit('ACTG');
is $a->repeat_unit, 'ACTG';

$a->repeat_count(3);
is $a->repeat_count, 3;

