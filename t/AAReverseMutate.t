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

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$obj = new Bio::Variation::AAReverseMutate
    ('-aa_ori' => 'F', 
     '-aa_mut' => 'S'
     );
test 2, defined $obj && ref($obj) =~ /Bio::Variation::AAReverseMutate/;


test 3, ($obj->aa_ori eq 'F' );

test 4, ($obj->aa_mut eq 'S' );

@points = $obj->each_Variant;
# F>S has two solutions
test 5, (scalar @points  == 2 );

$obj->codon_ori('ttc');
test 6, defined $obj;

#now there should be only one left
@points = $obj->each_Variant;
test 7, (scalar @points  == 1 );

$obj->codon_table(3);
test 8, ( $obj->codon_table == 3);

#Check the returned object
$rna = pop @points;
test 9, ( $rna->isa('Bio::Variation::RNAChange') );

test 10, ($rna->length == 1 );
test 11, ($rna->allele_ori->seq eq 't' );
test 12, ($rna->allele_mut->seq eq 'c' );


test 13, ($rna->codon_ori eq 'ttc' ), "Codon_ori is |". $rna->codon_ori. "|";

test 14, ($rna->codon_pos == 2 );

$obj->codon_table(11);
test 15, ( $obj->codon_table == 11);


