## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
# Test script for Bio::UnivAln.pm
# Georg Fuellen, adapted for the Bioperl distribution by Steve A. Chervitz
# Very rudimentary. Eventually will incorporate Georg's univaln.t2
#
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

BEGIN {	
    $| = 1; print "1..2\n"; 
}
END {
#   print "not ok 1\n" unless $loaded;
#   unlink $testout;  # commented out since you may want to check it...
}

use lib '../';
use Bio::UnivAln;

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}
my($s,@s);

test 1, $aln = Bio::UnivAln->new(-seqs=>"TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTG--TCAACTG");
test 2, print $aln->layout("fasta");

# print STDERR "\n\n\nTHIS TEST SCRIPT PREFORMS ONLY ONE BASIC TEST.\n";
# print STDERR "For intensive testing, see t/univaln.t2. Run that script via\n";
# print STDERR "% perl t/univaln.t2 > t/my_univaln\n";
# print STDERR "and compare the output with t/univaln.o\n";
# print STDERR "Expected error messages can be found in t/univaln2_expected_errors\n\n";
 
