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


my $seq = Bio::PrimarySeq->new(-seq=>'TTGGTGGCGTCAACT',
			-display_id => 'new-id',
			-moltype => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
print "ok 2\n";  

$seq->accession_number();
$seq->seq();
$seq->display_id();

print "ok 3\n";

$trunc = $seq->trunc(1,4);

print "ok 4\n";

if( $trunc->seq() ne 'TTGG' ) {
   print "not ok 5\n";
   $s = $trunc->seq();
   print STDERR "Expecting TTGG. Got $s\n";
} else {
   print "ok 5\n";
}

$rev = $seq->revcom();

print "ok 6\n";

if( $rev->seq() ne 'AGTTGACGCCACCAA' ) {
   print "not ok 7\n";
} else {
   print "ok 7\n";
}

#
# Translate
#

$aa = $seq->translate();

#print $seq->seq, "ok 8\n"; # you cannot simply comment out tests -- they
                            # will be treated as failed!

if( $aa->seq ne 'LVAST' ) {
    print STDERR "Translation: ", $aa->seq, "\n";
    print "not ok 8\n";
} else {
    print "ok 8\n";
}

$seq->seq('TTGGTGGCGTCAACTTAA');

#                    $stop, $unknown, $frame, $tableid, $fullCDS, $throw
$aa = $seq->translate(undef, undef, undef, undef, 1);

# tests for non-Methionin initiator codon (AGT) coding for M
if( $aa->seq ne 'MVAST' ) {
    print STDERR "Translation: ", $aa->seq, "\n";
    print "not ok 9\n";
} else {
    print "ok 9\n";
}

