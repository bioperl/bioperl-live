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
BEGIN { $| = 1; print "1..13\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::PrimarySeq;

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


my $seq = Bio::PrimarySeq->new(-seq=>'TTGGTGGCGTCAACT',
			-display_id => 'new-id',
			-moltype => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
test 2, defined $seq;

test 3, $seq->accession_number() eq 'X677667';
test 4, $seq->seq() eq 'TTGGTGGCGTCAACT';
test 5, $seq->display_id(), 'new-id';

test 6, $seq->display_id() eq 'new-id';
test 7, $seq->moltype() eq 'dna';

$trunc = $seq->trunc(1,4);
test 8, defined $trunc && $trunc->length == 4;
test 9, ( $trunc->seq() eq 'TTGG' ), "Expecting TTGG. Got ".$trunc->seq();

$rev = $seq->revcom();
test 10, defined $rev; 

test 11, ( $rev->seq() eq 'AGTTGACGCCACCAA' );

#
# Translate
#

$aa = $seq->translate();

test 12, ( $aa->seq eq 'LVAST' ), "Translation: ". $aa->seq;

$seq->seq('TTGGTGGCGTCAACTTAA');

$aa = $seq->translate(undef, undef, undef, undef, 1);

# tests for non-Methionin initiator codon (AGT) coding for M
test 13, ( $aa->seq eq 'MVAST' ), "Translation: ". $aa->seq;
