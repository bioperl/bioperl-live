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
BEGIN { $| = 1; print "1..15\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

#use lib '../';
use Bio::Variation::AAReverseMutate;
$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

$obj = new Bio::Variation::AAReverseMutate
    ('-aa_ori' => 'F', 
     '-aa_mut' => 'S'
     );
print "ok 2\n";  


if ($obj->aa_ori eq 'F' ) {
    print "ok 3\n";  
} else {
    print "not ok 3\n";
} 

if ($obj->aa_mut eq 'S' ) {
    print "ok 4\n";  
} else {
    print "not ok 4\n";
}

@points = $obj->each_Variant;
# F>S has two solutions
if (scalar @points  == 2 ) {
    print "ok 5\n";
} else {
    print "not ok 5\n";
} 

$obj->codon_ori('ttc');
print "ok 6\n";

#now there should be only one left
@points = $obj->each_Variant;
if (scalar @points  == 1 ) {
    print "ok 7\n";  
} else {
    print "not ok 7\n";
}

$obj->codon_table(3);
if( $obj->codon_table == 3) {
    print "ok 8\n";  
} else {
    print "not ok 8\n";
} 

#Check the returned object
$rna = pop @points;
if( $rna->isa('Bio::Variation::RNAChange') ) {
    print "ok 9\n";
} else {
    print "not ok 9\n";
} 

if ($rna->length == 1 ) {
    print "ok 10\n";  
} else {
    print "not ok 10\n";
} 

if ($rna->allele_ori->seq eq 't' ) {
    print "ok 11\n";  
} else {
    print "not ok 11\n";
} 

if ($rna->allele_mut->seq eq 'c' ) {
    print "ok 12\n";  
} else {
    print "not ok 12\n";
} 


if ($rna->codon_ori eq 'ttc' ) {
    print "ok 13\n";  
} else {
    print "Codon_ori is |", $rna->codon_ori, "|\n";
    print "not ok 13\n";
}

if ($rna->codon_pos == 2 ) {
    print "ok 14\n";  
} else {
    print "not ok 14\n";
} 

$obj->codon_table(11);
if( $obj->codon_table == 11) {
    print "ok 15\n";  
} else {
    print "not ok 15\n";
} 


