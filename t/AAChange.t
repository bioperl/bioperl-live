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


sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$obj = Bio::Variation::AAChange -> new;

test 2, defined $obj && ref($obj) =~ /Bio::Variation::AAChange/;

$obj->start(3);           
test 3, ($obj->start == 3 );


$obj->end(3); 
test 4, ($obj->end == 3 );

$obj->length(3);

test 5, ($obj->length == 3 );

$obj->strand('1');  
test 6, ($obj->strand eq '1' );

test 7, ($obj->primary_tag eq 'Variation' );

$obj->source_tag('source');
test 8, ($obj->source_tag eq 'source' );

$obj->frame(2);   
test 9, ($obj->frame ==2 );

$obj->score(2);   
test 10, ($obj->score ==2 );

$obj->isMutation(1); 
test 11, ($obj->isMutation );

$a1 = Bio::Variation::Allele->new(-seq => 'V');
$obj->allele_ori($a1);

test 12, ($obj->allele_ori->seq eq 'V' );

$a2 = Bio::Variation::Allele->new('-seq' => 'A');
$obj->add_Allele($a2);

test 13, ($obj->allele_mut->seq eq 'A' );

$obj->upStreamSeq('upStreamSeq'); 
test 14, ($obj->upStreamSeq eq 'upStreamSeq' );

$obj->dnStreamSeq('dnStreamSeq'); 
test 15, ($obj->dnStreamSeq eq 'dnStreamSeq' );

test 16, ($obj->label eq 'substitution' );

$obj->status('proven'); 
test 17, ($obj->status eq 'proven' );

$obj->proof('experimental'); 
test 18, ($obj->proof eq 'experimental' );

$obj->region('region'); 
test 19, ($obj->region eq 'region' );

$obj->region_value('region_value'); 
test 20, ($obj->region_value eq 'region_value' );

$obj->numbering('coding'); 
test 21,  ($obj->numbering eq 'coding' );

$obj2 = Bio::Variation::RNAChange -> new(-start => 7, 
					  -end => 7,
					  -cds_end => 100,
					  -codon_pos => 1,
					  -upStreamSeq => 'acgcgcgcgc',
					  -dnStreamSeq => 'acgcgcgcgc'
					  );
$obj2->label('missense');
$obj->RNAChange($obj2);

test 22, ( $obj->trivname eq 'V3A' ), "Trivial name is !". $obj->trivname. "!\n";

$obj->mut_number(2);
test 23, ( $obj->mut_number == 2 );
