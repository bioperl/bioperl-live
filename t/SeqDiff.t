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

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$obj = Bio::Variation::SeqDiff -> new;

test 2, 1;

$obj->id('id');           
test 3, ($obj->id eq 'id' );


$obj->sysname('sysname'); 
test 4, ($obj->sysname eq 'sysname' );

$obj->trivname('trivname');

test 5, ($obj->trivname eq 'trivname' );

$obj->chromosome('chr');  

test 6, ($obj->chromosome eq 'chr' );

$obj->description('desc');
test 7, ($obj->description eq 'desc' );

$obj->numbering('numbering');
test 8, ($obj->numbering eq 'numbering' );

$obj->offset(100);   
test 9, ($obj->offset == 100 );

$obj->dna_ori('dna_ori'); 
test 10, ($obj->dna_ori eq 'dna_ori' );

$obj->dna_mut('dna_mut'); 
test 11, ($obj->dna_mut eq 'dna_mut' );

$obj->rna_ori('rna_ori'); 
test 12, ($obj->rna_ori eq 'rna_ori' );

$obj->rna_mut('rna_mut'); 
test 13, ($obj->rna_mut eq 'rna_mut' );

$obj->aa_ori('aa_ori'); 
test 14, ($obj->aa_ori eq 'aa_ori' );

$obj->aa_mut('aa_mut'); 
test 15, ($obj->aa_mut eq 'aa_mut' );

$m = Bio::Variation::DNAMutation->new;

test 16, 1;

$obj->add_Variant($m);
test 17, 1;

foreach $mm ( $obj->each_Variant ) {
    $mm->primary_tag('a');    
}
test 18, 1;

$obj->gene_symbol('fos');
test 19, ($obj->gene_symbol eq 'fos' );
