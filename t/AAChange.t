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
BEGIN { $| = 1; print "1..23\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Variation::Allele;
use Bio::Variation::AAChange;
use Bio::Variation::RNAChange;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


$obj = Bio::Variation::AAChange -> new;

print "ok 2\n";  

$obj->start(3);           
if ($obj->start == 3 ) {
    print "ok 3\n";  
} else {
    print "not ok 3\n";
} 


$obj->end(3); 
if ($obj->end == 3 ) {
    print "ok 4\n";  
} else {
    print "not ok 4\n";
} 

$obj->length(3);

if ($obj->length == 3 ) {
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

$obj->isMutation(1); 
if ($obj->isMutation ) {
    print "ok 11\n";  
} else {
    print "not ok 11\n";
}

$a1 = Bio::Variation::Allele->new(-seq => 'V');
$obj->allele_ori($a1);

if ($obj->allele_ori->seq eq 'V' ) {
    print "ok 12\n";  
} else {
    print "not ok 12\n";
} 


$a2 = Bio::Variation::Allele->new('-seq' => 'A');
$obj->add_Allele($a2);

if ($obj->allele_mut->seq eq 'A' ) {
    print "ok 13\n";  
} else {
    print "not ok 13\n";
}

$obj->upStreamSeq('upStreamSeq'); 
if ($obj->upStreamSeq eq 'upStreamSeq' ) {
    print "ok 14\n";  
} else {
    print "not ok 14\n";
} 

$obj->dnStreamSeq('dnStreamSeq'); 
if ($obj->dnStreamSeq eq 'dnStreamSeq' ) {
    print "ok 15\n";  
} else {
    print "not ok 15\n";
} 


if ($obj->label eq 'substitution' ) {
    print "ok 16\n";  
} else {
    print "not ok 16\n";
} 


$obj->status('proven'); 
if ($obj->status eq 'proven' ) {
    print "ok 17\n";  
} else {
    print "not ok 17\n";
} 


$obj->proof('experimental'); 
if ($obj->proof eq 'experimental' ) {
    print "ok 18\n";  
} else {
    print "not ok 18\n";
} 


$obj->region('region'); 
if ($obj->region eq 'region' ) {
    print "ok 19\n";  
} else {
    print "not ok 19\n";
} 


$obj->region_value('region_value'); 
if ($obj->region_value eq 'region_value' ) {
    print "ok 20\n";  
} else {
    print "not ok 20\n";
} 


$obj->numbering('coding'); 
if ($obj->numbering eq 'coding' ) {
    print "ok 21\n";  
} else {
    print "not ok 21\n";
} 


$obj2 = Bio::Variation::RNAChange -> new(-start => 7, 
					  -end => 7,
					  -cds_end => 100,
					  -codon_pos => 1,
					  -upStreamSeq => 'acgcgcgcgc',
					  -dnStreamSeq => 'acgcgcgcgc'
					  );
$obj2->label('missense');
$obj->RNAChange($obj2);

if ( $obj->trivname eq 'V3A' ) {
    print "ok 22\n";  
} else {
    print "Trivial name is !", $obj->trivname, "!\n";
    print "not ok 22\n";
} 

$obj->mut_number(2);
if ( $obj->mut_number == 2 ) {
    print "ok 23\n";  
} else {
    print "not ok 23\n";
} 
