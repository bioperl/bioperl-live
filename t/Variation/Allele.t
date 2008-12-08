# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
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
is $trunc->seq(), 'ACTG';

ok defined($rev = $a->revcom());
is $rev->seq(), 'CAGTCAGTCAGT';

$a->is_reference(1);
ok $a->is_reference;

$a->repeat_unit('ACTG');
is $a->repeat_unit, 'ACTG';

$a->repeat_count(3);
is $a->repeat_count, 3;
