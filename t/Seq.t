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
use Bio::Seq;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
			-moltype => 'dna' );
print "ok 2\n"; 

$trunc = $seq->trunc(1,4);

print "ok 3\n";

if( $trunc->seq() ne 'ACTG' ) {
   print "not ok 4\n";
} else {
   print "ok 4\n";
}

$trans = $seq->translate();

if( $trans->seq() ne 'TVAST' ) {
   print "not ok 5\n";
} else {
   print "ok 5\n";
}

# test ability to get str function

$t = $seq->seq();
if( $t eq 'ACTGTGGCGTCAACT' ) {
  print "ok 6\n";
}

$seq = Bio::Seq->new(-seq=>'actgtggcgtcaact',
		     -desc=>'Sample Bio::Seq object',
		     -display_id => 'something',
		     -accession_number => 'accnum',
		     -moltype => 'dna' );
print "ok 7\n"; 


$trans = $seq->translate();

if( $trans->seq() ne 'TVAST' ) {
   print "not ok 8\n";
} else {
   print "ok 8\n";
}

# basic methods

if( $seq->id() ne 'something' || $seq->accession_number ne 'accnum' ) {
    print "not ok 9\n";
    print "saw ",$seq->id,":",$seq->accession_number,":",$seq->primary_id,"\n";
} else {
  print "ok 9\n";
}

