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
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 37;
	use_ok('Bio::Variation::DNAMutation');
	use_ok('Bio::Variation::Allele');
}


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my($obj,$a1,$a2,$obj2);
$obj = Bio::Variation::DNAMutation -> new;

ok defined $obj;

$obj->start(3);           
is $obj->start, 3;

$obj->end(3); 
is $obj->end, 3;

$obj->length(2);
is $obj->length, 2;

$obj->strand('1');  
is $obj->strand, '1';

is $obj->primary_tag, 'Variation';

$obj->source_tag('source');
is $obj->source_tag, 'source';

$obj->frame(2);   
is $obj->frame,2;

$obj->score(2);   
is $obj->score, 2;

if( $obj->can('dna_mut') ) {
#test gff string
    $obj->dna_mut('dna_mut'); 
    is( $obj->dna_mut,'dna_mut');
}

$a1 = Bio::Variation::Allele->new(-seq => 'c');
$obj->allele_ori($a1);
 
is $obj->allele_ori->seq, 'c';

$a2 = Bio::Variation::Allele->new('-seq' => 'g');
$obj->allele_mut($a2);

is $obj->allele_mut->seq, 'g';

$obj->upStreamSeq('agcacctcccggcgccagtttgctg'); 
is $obj->upStreamSeq, 'agcacctcccggcgccagtttgctg';

$obj->dnStreamSeq('tgctgcagcagcagcagcagcagca'); 
is $obj->dnStreamSeq, 'tgctgcagcagcagcagcagcagca';


is $obj->label, 'point, transversion' ;

$obj->status('proven'); 
is $obj->status, 'proven';


$obj->proof('experimental'); 
is $obj->proof, 'experimental';


is $obj->restriction_changes, '-BbvI, +BstXI, -Fnu4HI, -TseI';

$obj->region('region'); 
is $obj->region, 'region';

$obj->region_value('region_value'); 
is $obj->region_value, 'region_value';

$obj->region_dist(-5); 
is $obj->region_dist, -5;

$obj->numbering('coding'); 
is $obj->numbering, 'coding';

ok not $obj->CpG;

$obj->mut_number(2);
is $obj->mut_number, 2;


ok defined ($obj2 = Bio::Variation::DNAMutation -> new
	    ('-mut_number' => 2));

is $obj2->mut_number, 2;


$obj->isMutation(1); 
ok $obj->isMutation;

$obj->add_Allele($a1);
$obj->add_Allele($a2);

is scalar ($obj->each_Allele), 2;


$obj = Bio::Variation::DNAMutation->new
    ('-start'         => 23,
     '-end'           => 24,
     '-length'        => 2,
     '-upStreamSeq'   => 'gt',
     '-dnStreamSeq'   => 'at',
     '-proof'         => 'experimental',
     '-isMutation'    => 1,
     '-mut_number'    => 2
     );

is $obj->start(), 23;
is $obj->end(), 24;
is $obj->length(), 2;
is $obj->upStreamSeq(), 'gt';
is $obj->dnStreamSeq(), 'at';
is $obj->proof(), 'experimental';
is $obj->mut_number(), 2;
ok $obj->isMutation;
