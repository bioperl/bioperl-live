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
BEGIN { $| = 1; print "1..3\n"; }
END {print "not ok 1\n" unless $loaded;}

use Bio::AnnSeqIO;
use Bio::AnnSeq;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

$ast = Bio::AnnSeqIO->new( '-format' => 'EMBL' , -file => 't/roa1.dat');

while( my $as = $ast->next_annseq() ) {
       if( ! defined $as->seq ) {
	   print "not ok 2\n";
	   }
      }


print "ok 2\n";

$ast = Bio::AnnSeqIO->new( '-format' => 'GenBank' , -file => 't/roa1.genbank');

while( my $as = $ast->next_annseq() ) {
       if( ! defined $as->seq ) {
	   print "not ok 3\n";
	   }
      }


print "ok 3\n";
