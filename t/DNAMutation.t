# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..28\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Variation::DNAMutation;
use Bio::Variation::Allele;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$obj = Bio::Variation::DNAMutation -> new;

test 2, defined $obj;

$obj->start(3);           
test 3, ($obj->start == 3 );

$obj->end(3); 
test 4, ($obj->end == 3 );

$obj->length(2);

test 5, ($obj->length == 2 );

$obj->strand('1');  
test 6, ($obj->strand eq '1' );

test 7, ($obj->primary_tag eq 'Variation' );

$obj->source_tag('source');
test 8, ($obj->source_tag eq 'source' );

$obj->frame(2);   
test 9,  ($obj->frame ==2 );

$obj->score(2);   
test 10, ($obj->score ==2 );

#test gff string
#$obj->dna_mut('dna_mut'); 
#if ($obj->dna_mut eq 'dna_mut' ) {
#    print "ok 11\n";  
#} else {
#    print "not ok 11\n";
#} 
test 11, 1;

$a1 = Bio::Variation::Allele->new(-seq => 'c');
$obj->allele_ori($a1);
 
test 12, ($obj->allele_ori->seq eq 'c' );

$a2 = Bio::Variation::Allele->new('-seq' => 'g');
$obj->allele_mut($a2);

test 13, ($obj->allele_mut->seq eq 'g' );

$obj->upStreamSeq('agcacctcccggcgccagtttgctg'); 
test 14, ($obj->upStreamSeq eq 'agcacctcccggcgccagtttgctg' );

$obj->dnStreamSeq('tgctgcagcagcagcagcagcagca'); 
test 15, ($obj->dnStreamSeq eq 'tgctgcagcagcagcagcagcagca' );


test 16, ($obj->label eq 'point, transversion' );

$obj->status('proven'); 
test 17, ($obj->status eq 'proven' );


$obj->proof('experimental'); 
test 18, ($obj->proof eq 'experimental' );


test 19, ($obj->restriction_changes eq '-BbvI, +BstXI, -Fnu4HI, -TseI' );

$obj->region('region'); 
test 20, ($obj->region eq 'region' );

$obj->region_value('region_value'); 
test 21, ($obj->region_value eq 'region_value' );

$obj->numbering('coding'); 
test 22, ($obj->numbering eq 'coding' );

test 23, (not $obj->CpG );

$obj->mut_number(2);
test 24, ( $obj->mut_number == 2 );


$obj2 = Bio::Variation::DNAMutation -> new
	 ('-mut_number' => 2);

test 25, ( $obj2->mut_number == 2 );


$obj->isMutation(1); 
test 26, ($obj->isMutation );

$obj->add_Allele($a1);
test 27, defined $obj;

$obj->add_Allele($a2);
test 28, (scalar ($obj->each_Allele) == 2 );
