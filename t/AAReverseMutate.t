# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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

    plan tests => 16;
}

use Bio::Variation::AAReverseMutate;
ok(1);

my $obj = new Bio::Variation::AAReverseMutate
    ('-aa_ori' => 'F', 
     '-aa_mut' => 'S'
     );
ok defined $obj;
ok $obj->isa('Bio::Variation::AAReverseMutate');

ok $obj->aa_ori, 'F';

ok $obj->aa_mut, 'S';

my @points = $obj->each_Variant;
# F>S has two solutions
ok scalar @points, 2;

$obj->codon_ori('ttc');
ok defined $obj;

#now there should be only one left
@points = $obj->each_Variant;
ok scalar @points, 1;

$obj->codon_table(3);
ok $obj->codon_table, 3;

#Check the returned object
my $rna = pop @points;
ok $rna->isa('Bio::Variation::RNAChange');

ok $rna->length, 1;
ok $rna->allele_ori->seq, 't';
ok $rna->allele_mut->seq, 'c';


ok $rna->codon_ori, 'ttc', "Codon_ori is |". $rna->codon_ori. "|";

ok $rna->codon_pos, 2;

$obj->codon_table(11);
ok $obj->codon_table, 11;
