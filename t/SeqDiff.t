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
BEGIN { $| = 1; print "1..19\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


$obj = Bio::Variation::SeqDiff -> new;

print "ok 2\n";  

$obj->id('id');           
if ($obj->id eq 'id' ) {
    print "ok 3\n";  
} else {
    print "not ok 3\n";
} 


$obj->sysname('sysname'); 
if ($obj->sysname eq 'sysname' ) {
    print "ok 4\n";  
} else {
    print "not ok 4\n";
} 

$obj->trivname('trivname');

if ($obj->trivname eq 'trivname' ) {
    print "ok 5\n";  
} else {
    print "not ok 5\n";
} 

$obj->chromosome('chr');  

if ($obj->chromosome eq 'chr' ) {
    print "ok 6\n";  
} else {
    print "not ok 6\n";
} 

$obj->description('desc');
if ($obj->description eq 'desc' ) {
    print "ok 7\n";  
} else {
    print "not ok 7\n";
} 

$obj->numbering('numbering');
if ($obj->numbering eq 'numbering' ) {
    print "ok 8\n";  
} else {
    print "not ok 8\n";
} 

$obj->offset(100);   
if ($obj->offset == 100 ) {
    print "ok 9\n";  
} else {
    print "not ok 9\n";
} 

$obj->dna_ori('dna_ori'); 
if ($obj->dna_ori eq 'dna_ori' ) {
    print "ok 10\n";  
} else {
    print "not ok 10\n";
} 

$obj->dna_mut('dna_mut'); 
if ($obj->dna_mut eq 'dna_mut' ) {
    print "ok 11\n";  
} else {
    print "not ok 11\n";
} 

$obj->rna_ori('rna_ori'); 
if ($obj->rna_ori eq 'rna_ori' ) {
    print "ok 12\n";  
} else {
    print "not ok 12\n";
} 

$obj->rna_mut('rna_mut'); 
if ($obj->rna_mut eq 'rna_mut' ) {
    print "ok 13\n";  
} else {
    print "not ok 13\n";
}

$obj->aa_ori('aa_ori'); 
if ($obj->aa_ori eq 'aa_ori' ) {
    print "ok 14\n";  
} else {
    print "not ok 14\n";
} 

$obj->aa_mut('aa_mut'); 
if ($obj->aa_mut eq 'aa_mut' ) {
    print "ok 15\n";  
} else {
    print "not ok 15\n";
} 

$m = Bio::Variation::DNAMutation->new;

print "ok 16\n";  

$obj->add_Variant($m);
print "ok 17\n";  

foreach $mm ( $obj->each_Variant ) {
    $mm->primary_tag('a');
    
}
print  "ok 18\n";       

$obj->gene_symbol('fos');
if ($obj->gene_symbol eq 'fos' ) {
    print "ok 19\n";  
} else {
    print "not ok 19\n";
} 
