## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
# Test script for Bio::Tools::Fasta.pm
# Steve A. Chervitz, sac@neomorphic.com
# Fairly rudimentary.
# Strategy for this test was borrowed from L. Stein's BoulderIO package.

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
BEGIN { $| = 1; print "1..5\n"; }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::Fasta;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my @seqs;

test 2, $fasta = new Bio::Tools::Fasta (-file => 't/seqs.fas',
	                                -seqs => 1,
		                        -save_array => \@seqs,
		                        -parse => 1,
	                                -edit_id => 1,
	                               );

test 3, scalar(@seqs) == 6, "Number of seqs = ${\scalar(@seqs)}";

test 4, $fasta->num_seqs == 6, "Number of seqs = ${\$fasta->num_seqs}";

print "First sequence:\n";

test 5, print $seqs[0]->layout('fasta');





