# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;

BEGIN { plan tests => 27 }

use Bio::Variation::DNAMutation;
use Bio::Variation::Allele;


ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my($obj,$a1,$a2,$obj2);
$obj = Bio::Variation::DNAMutation -> new;

ok defined $obj;

$obj->start(3);           
ok $obj->start, 3;

$obj->end(3); 
ok $obj->end, 3;

$obj->length(2);
ok $obj->length, 2;

$obj->strand('1');  
ok $obj->strand, '1';

ok $obj->primary_tag, 'Variation';

$obj->source_tag('source');
ok $obj->source_tag, 'source';

$obj->frame(2);   
ok $obj->frame,2;

$obj->score(2);   
ok $obj->score, 2;

if( $obj->can('dna_mut') ) {
#test gff string
    $obj->dna_mut('dna_mut'); 
    ok( $obj->dna_mut,'dna_mut');
}

$a1 = Bio::Variation::Allele->new(-seq => 'c');
$obj->allele_ori($a1);
 
ok $obj->allele_ori->seq, 'c';

$a2 = Bio::Variation::Allele->new('-seq' => 'g');
$obj->allele_mut($a2);

ok $obj->allele_mut->seq, 'g';

$obj->upStreamSeq('agcacctcccggcgccagtttgctg'); 
ok $obj->upStreamSeq, 'agcacctcccggcgccagtttgctg';

$obj->dnStreamSeq('tgctgcagcagcagcagcagcagca'); 
ok $obj->dnStreamSeq, 'tgctgcagcagcagcagcagcagca';


ok $obj->label, 'point, transversion' ;

$obj->status('proven'); 
ok $obj->status, 'proven';


$obj->proof('experimental'); 
ok $obj->proof, 'experimental';


ok $obj->restriction_changes, '-BbvI, +BstXI, -Fnu4HI, -TseI';

$obj->region('region'); 
ok $obj->region, 'region';

$obj->region_value('region_value'); 
ok $obj->region_value, 'region_value';

$obj->numbering('coding'); 
ok $obj->numbering, 'coding';

ok not $obj->CpG;

$obj->mut_number(2);
ok $obj->mut_number, 2;


ok defined ($obj2 = Bio::Variation::DNAMutation -> new
	    ('-mut_number' => 2));

ok $obj2->mut_number, 2;


$obj->isMutation(1); 
ok $obj->isMutation;

$obj->add_Allele($a1);
$obj->add_Allele($a2);

ok scalar ($obj->each_Allele), 2;
