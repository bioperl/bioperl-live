## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


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
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Index::Fasta;
use Bio::Index::SwissPfam;
use Bio::Index::EMBL;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

chomp( $dir = `pwd` );
$ind = Bio::Index::Fasta->new('Wibbl', 'WRITE');
$ind->make_index("$dir/t/seqs.fas");

print "ok 2\n";
$ind = 0;

$ind = Bio::Index::SwissPfam->new('Wibbl2', 'WRITE');
$ind->make_index("$dir/t/swisspfam.data");

print "ok 3\n";
$ind = 0;

# don't test EMBL yet. Bad...

system("rm -f Wibbl*");




