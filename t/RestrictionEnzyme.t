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
BEGIN { $| = 1; print "1..24\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::RestrictionEnzyme;
use Bio::Seq;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    if( !defined $msg ) {
	$msg = "No message";
	}
	 
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

# pre declare 6cutters to prevent warnings

@sixcutters = ();

$dna = 'CCTCCGGGGACTGCCGTGCCGGGCGGGAATTCGCCATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGA';

# Build sequence and restriction enzyme objects.
test 2, $seq = new Bio::Seq(-ID  =>'test_seq',
			     -SEQ =>$dna); 

test 3, $seq->id();
test 4, $re  = new Bio::Tools::RestrictionEnzyme(-NAME=>'EcoRI');
test 5, $re->seq->id;
test 6, $re->site;
test 7, $re->palindromic;
test 8, @fragments = $re->cut_seq($seq), "can't cut sequence";
test 9, scalar @fragments == 2;

test 10, ($l1 = length $fragments[0]) == 27;
test 11, ($l2 = length $fragments[1]) == 43;
test 12, $l1 + $l2 == $seq->length, "sum of frag lengths != length of original seq\n";

test 13, @sixcutters = $re->available_list(6), "can't get list of 6-cutters";
test 14, $re->is_available('HindIII');

test 15, @positions = $re->cuts_seq_at($seq);
test 16, scalar(@positions) == 1;
test 17, $positions[0] == 27;

# TspRI cuts at the pattern: CASTGNN, after the second N.
# forward strand cuts (after the bar)
#                   |                                               |
$dna = "ACGTGGCAGTGATGGCGGGAATTCGCCATGGCGACCCTGGAAAAGCTGATGATTCAGTGAAAACTG";
#          |                                               |
# reverse strand cuts (after the bar)
# Build sequence and restriction enzyme objects.

$seq = new Bio::Seq(-ID  =>'test_seq', -SEQ =>$dna); 
	
test 18, $re  = new Bio::Tools::RestrictionEnzyme(-NAME=>'TspRI');
test 19, @positions = $re->cuts_seq_at($seq);
test 20, scalar(@positions) == 4;
test 21, $positions[0] == 4;	# this one's on the reverse strand.
test 22, $positions[1] == 13;	# this one's on the forward strand.
test 23, $positions[2] == 52;	# this one's on the reverse strand.
test 24, $positions[3] == 61;	# this one's on the forward strand.



