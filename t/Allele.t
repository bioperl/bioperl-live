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
use Bio::Variation::Allele;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


my $a = Bio::Variation::Allele->new(-seq=>'ACTGACTGACTG',
			-display_id => 'new-id',
			-moltype => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
print "ok 2\n";  

$a->accession_number();
$a->seq();
$a->display_id();

print "ok 3\n";

$trunc = $a->trunc(1,4);

print "ok 4\n";

if( $trunc->seq() ne 'ACTG' ) {
   print "not ok 5\n";
   $s = $trunc->seq();
   print STDERR "Expecting ACTG. Got $s\n";
} else {
   print "ok 5\n";
}

$rev = $a->revcom();

print "ok 6\n";

if( $rev->seq() ne 'CAGTCAGTCAGT' ) {
   print "not ok 7\n";
} else {
   print "ok 7\n";
}

$a->is_reference(1);
print "ok 8\n";

if( $a->is_reference ) {
    print "ok 9\n";
} else {
    print "not ok 9\n";
}

$a->repeat_unit('ACTG');
print "ok 10\n";

if( $a->repeat_unit eq 'ACTG' ) {
    print "ok 11\n";
} else {
    print "not ok 11\n";
}


$a->repeat_count(3);
print "ok 12\n";

if( $a->repeat_count == 3 ) {
   print "ok 13\n";
} else {
   print "not ok 13\n";
}

