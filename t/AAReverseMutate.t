# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
## We start with some black magic to print on failure.

use Test;
use strict;
BEGIN { plan tests => 16 }

use Bio::Variation::AAReverseMutate;
ok(1);
## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


my $obj = new Bio::Variation::AAReverseMutate
    ('-aa_ori' => 'F', 
     '-aa_mut' => 'S'
     );
ok defined $obj;
ok ref($obj), qr/Bio::Variation::AAReverseMutate/;

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


