## Bioperl Test Harness Script for Modules
##


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
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Tools::HMMER::Domain;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

$domain = Bio::Tools::HMMER::Domain->new();

$domain->start(50);
$domain->end(200);
$domain->start_hmm(10);
$domain->end_hmm(100);
$domain->seqbits(50);
$domain->bits(20);
$domain->evalue(0.0001);
$domain->seqname('silly');

print "ok 2\n";

# test that we can get out forward and reverse homol_SeqFeatures
$homol = $domain->homol_SeqFeature();
if( $homol->start() == 10 ) {
    print "ok 3\n";
}

$rev = $homol->homol_SeqFeature();
if( $rev->start() == 50 ) {
    print "ok 4\n"
}

