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

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}
 
$obj = Bio::Variation::RNAChange -> new;

test 2, defined $obj;

$obj->start(4);           
test 3, ($obj->start == 4 );

$obj->end(4); 
test 4, ($obj->end == 4 );

$obj->length(1);

test 5, ($obj->length == 1 );

$obj->strand('1');  
test 6, ($obj->strand eq '1' );

test 7, ($obj->primary_tag eq 'Variation' );

$obj->source_tag('source');
test 8, ($obj->source_tag eq 'source' );

$obj->frame(2);   
test 9, ($obj->frame ==2 );

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

$a1 = Bio::Variation::Allele->new(-seq => 'g');
$obj->allele_ori($a1);

test 12, ($obj->allele_ori->seq eq 'g' );

$a2 = Bio::Variation::Allele->new('-seq' => 'a');
$obj->allele_mut($a2);

test 13, ($obj->allele_mut->seq eq 'a' );

$obj->upStreamSeq('gaagattcagccaagctcaaggatg'); 
test 14, ($obj->upStreamSeq eq 'gaagattcagccaagctcaaggatg' );

$obj->cds_end(1000); 
test 15, ($obj->cds_end == 1000 );

$obj->dnStreamSeq('aagtgcagttagggctgggaagggt'); 
test 16, ($obj->dnStreamSeq eq 'aagtgcagttagggctgggaagggt' );

$obj->codon_pos(1); 
test 17, ($obj->codon_pos == 1 );

$obj3 = Bio::Variation::AAChange -> new;
$obj3->start(2);
$obj->AAChange($obj3);

test 18, ($obj->label eq 'missense' ), "label is". $obj->label;


$obj->status('proven'); 
test 19, ($obj->status eq 'proven' );

$obj->proof('experimental'); 
test 20, ($obj->proof eq 'experimental' );
test 21, ($obj->restriction_changes eq '-BccI' );

$obj->region('coding'); 
test 22, ($obj->region eq 'coding' );
$obj->numbering('coding'); 
test 23, ($obj->numbering eq 'coding' );

test 24, ($obj->codon_ori eq 'gaa' ), "Codon_ori is |". $obj->codon_ori. "|";

test 25, ($obj->codon_mut eq 'aaa' ), "Codon_mut is |". $obj->codon_mut. "|";


$obj->codon_pos(1); 
test 26, ($obj->codon_pos == 1 );
test 27, ( $obj->codon_table == 1 );

$obj->codon_table(3);
test 28, ( $obj->codon_table == 3 );

$obj->mut_number(2);
test 29, ( $obj->mut_number == 2 );

$obj->verbose(2);
test 30, ( $obj->verbose == 2 );


