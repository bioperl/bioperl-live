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
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Tools::HMMER::Domain;
use Bio::Tools::HMMER::Set;
use Bio::Tools::HMMER::Results;

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
$homol = $domain->feature2();
if( $homol->start() == 10 ) {
    print "ok 3\n";
} else {
   print "not ok 3\n";
}


$rev = $domain;

if( $rev->start() == 50 ) {
    print "ok 4\n"
}else {
   print "not ok 4\n";
}


$set = Bio::Tools::HMMER::Set->new();
$set->add_Domain($domain);

@doms = $set->each_Domain();
$dom = shift @doms;
if( $dom->start() == 50 ) {
   print "ok 5\n";
} else {
   print "not ok 5\n";
}

$set->bits(300);
$set->evalue(0.0001);
$set->name('sillyname');

print "ok 6\n";

$res = Bio::Tools::HMMER::Results->new( -file => 't/hmmsearch.out' , -type => 'hmmsearch');
foreach $set ( $res->each_Set) {
  foreach $domain ( $set->each_Domain ) {
    #print STDERR "Got domain ",$domain->seqname," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
  }
}

print "ok 7\n";

if( $res->number != 1215 ) {
   print STDERR "\nBad number of domains. Expecting 1215. Got",$res->number,"\n";
   print "not ok 8\n";
} else {
   print "ok 8\n";
}

$res = Bio::Tools::HMMER::Results->new( -file => 't/hmmpfam.out' , -type => 'hmmpfam');
if( $res->number != 2 ) {
   print STDERR "\nBad number of domains. Expecting 2. Got",$res->number,"\n";
   print "not ok 9\n";
} else {
   print "ok 9\n";
}








