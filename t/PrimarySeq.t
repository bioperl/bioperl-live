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
BEGIN { $| = 1; print "1..9\n"; 
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

$seq->accession_number();
$seq->seq();
$seq->display_id();

test 3, $seq->display_id() eq 'new-id';

$trunc = $seq->trunc(1,4);
test 4, defined $trunc;
test 5, ( $trunc->seq() eq 'TTGG' ), "Expecting TTGG. Got ".$trunc->seq();

$rev = $seq->revcom();
test 6, defined $rev; 

test 7, ( $rev->seq() eq 'AGTTGACGCCACCAA' );

#
# Translate
#

$aa = $seq->translate();

#print $seq->seq, "ok 8\n"; # you cannot simply comment out tests -- they
                            # will be treated as failed!

test 8, ( $aa->seq eq 'LVAST' ), "Translation: ". $aa->seq;

$seq->seq('TTGGTGGCGTCAACTTAA');

#                    $stop, $unknown, $frame, $tableid, $fullCDS, $throw
$aa = $seq->translate(undef, undef, undef, undef, 1);

# tests for non-Methionin initiator codon (AGT) coding for M
test 9, ( $aa->seq eq 'MVAST' ), "Translation: ". $aa->seq;
