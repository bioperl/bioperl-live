# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 16);
	
	use_ok('Bio::Variation::AAReverseMutate');
}

my $obj = Bio::Variation::AAReverseMutate->new
    ('-aa_ori' => 'F', 
     '-aa_mut' => 'S'
     );
ok defined $obj;
isa_ok($obj, 'Bio::Variation::AAReverseMutate');

is $obj->aa_ori, 'F';

is $obj->aa_mut, 'S';

my @points = $obj->each_Variant;
# F>S has two solutions
is scalar @points, 2;

$obj->codon_ori('ttc');
ok defined $obj;

#now there should be only one left
@points = $obj->each_Variant;
is scalar @points, 1;

$obj->codon_table(3);
is $obj->codon_table, 3;

#Check the returned object
my $rna = pop @points;
isa_ok($rna, 'Bio::Variation::RNAChange');

is $rna->length, 1;
is $rna->allele_ori->seq, 't';
is $rna->allele_mut->seq, 'c';

is $rna->codon_ori, 'ttc', "Codon_ori is |". $rna->codon_ori. "|";

is $rna->codon_pos, 2;

$obj->codon_table(11);
is $obj->codon_table, 11;
