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
BEGIN { $| = 1; print "1..30\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Variation::Allele;
use Bio::Variation::RNAChange;
use Bio::Variation::DNAMutation;
use Bio::Variation::AAChange;
$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 
 
$obj = Bio::Variation::RNAChange -> new;

print "ok 2\n";  

$obj->start(4);           
if ($obj->start == 4 ) {
    print "ok 3\n";  
} else {
    print "not ok 3\n";
} 


$obj->end(4); 
if ($obj->end == 4 ) {
    print "ok 4\n";  
} else {
    print "not ok 4\n";
} 

$obj->length(1);

if ($obj->length == 1 ) {
    print "ok 5\n";  
} else {
    print "not ok 5\n";
} 

$obj->strand('1');  
if ($obj->strand eq '1' ) {
    print "ok 6\n";  
} else {
    print "not ok 6\n";
} 

if ($obj->primary_tag eq 'Variation' ) {
    print "ok 7\n";
} else {
    print "not ok 7\n";
} 

$obj->source_tag('source');
if ($obj->source_tag eq 'source' ) {
    print "ok 8\n";  
} else {
    print "not ok 8\n";
} 

$obj->frame(2);   
if ($obj->frame ==2 ) {
    print "ok 9\n";  
} else {
    print "not ok 9\n";
} 

$obj->score(2);   
if ($obj->score ==2 ) {
    print "ok 10\n";  
} else {
    print "not ok 10\n";
} 

#test gff string
#$obj->dna_mut('dna_mut'); 
#if ($obj->dna_mut eq 'dna_mut' ) {
#    print "ok 11\n";  
#} else {
#    print "not ok 11\n";
#} 
print "ok 11\n";  

$a1 = Bio::Variation::Allele->new(-seq => 'g');
$obj->allele_ori($a1);

if ($obj->allele_ori->seq eq 'g' ) {
    print "ok 12\n";  
} else {
    print "not ok 12\n";
} 

$a2 = Bio::Variation::Allele->new('-seq' => 'a');
$obj->allele_mut($a2);

if ($obj->allele_mut->seq eq 'a' ) {
    print "ok 13\n";  
} else {
    print "not ok 13\n";
}

$obj->upStreamSeq('gaagattcagccaagctcaaggatg'); 
if ($obj->upStreamSeq eq 'gaagattcagccaagctcaaggatg' ) {
    print "ok 14\n";  
} else {
    print "not ok 14\n";
} 


$obj->cds_end(1000); 
if ($obj->cds_end == 1000 ) {
    print "ok 15\n";  
} else {
    print "not ok 15\n";
} 


$obj->dnStreamSeq('aagtgcagttagggctgggaagggt'); 
if ($obj->dnStreamSeq eq 'aagtgcagttagggctgggaagggt' ) {
    print "ok 16\n";  
} else {
    print "not ok 16\n";
} 

$obj->codon_pos(1); 
if ($obj->codon_pos == 1 ) {
    print "ok 17\n";  
} else {
    print "not ok 17\n";
} 

$obj3 = Bio::Variation::AAChange -> new;
$obj3->start(2);
$obj->AAChange($obj3);

if ($obj->label eq 'missense' ) {
    print "ok 18\n";  
} else {
    print "label is", $obj->label, "\n";
    print "not ok 18\n" ;
} 


$obj->status('proven'); 
if ($obj->status eq 'proven' ) {
    print "ok 19\n";  
} else {
    print "not ok 19\n";
} 


$obj->proof('experimental'); 
if ($obj->proof eq 'experimental' ) {
    print "ok 20\n";  
} else {
    print "not ok 20\n";
} 


if ($obj->restriction_changes eq '-BccI' ) {
    print "ok 21\n";  
} else {
    print "not ok 21\n";
}


$obj->region('coding'); 
if ($obj->region eq 'coding' ) {
    print "ok 22\n";  
} else {
    print "not ok 22\n";
}

$obj->numbering('coding'); 
if ($obj->numbering eq 'coding' ) {
    print "ok 23\n";  
} else {
    print "not ok 23\n";
} 

if ($obj->codon_ori eq 'gaa' ) {
    print "ok 24\n";  
} else {
    print "Codon_ori is |", $obj->codon_ori, "|\n";
    print "not ok 24\n";
} 

if ($obj->codon_mut eq 'aaa' ) {
    print "ok 25\n";  
} else {
    print "Codon_mut is |", $obj->codon_mut, "|\n";
    print "not ok 25\n";
}


$obj->codon_pos(1); 
if ($obj->codon_pos == 1 ) {
    print "ok 26\n";  
} else {
    print "not ok 26\n";
} 

if ( $obj->codon_table == 1 ) {
    print "ok 27\n";  
} else {
    print "not ok 27\n";
} 

$obj->codon_table(3);
if ( $obj->codon_table == 3 ) {
    print "ok 28\n";  
} else {
    print "not ok 28\n";
} 

$obj->mut_number(2);
if ( $obj->mut_number == 2 ) {
    print "ok 29\n";  
} else {
    print "not ok 29\n";
} 

$obj->verbose(2);
if ( $obj->verbose ) {
    print "ok 30\n";  
} else {
    print "not ok 30\n";
} 


