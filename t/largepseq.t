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
BEGIN { $| = 1; print "1..3\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}
$loaded = 1;
print "ok 1\n";
use lib '../';

use Bio::Seq::LargePrimarySeq;


$pseq = Bio::Seq::LargePrimarySeq->new();

$pseq->add_sequence_as_string('ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAAT');
$pseq->add_sequence_as_string('GTTTGGGGTTAAACCCCTTTGGGGGGT');

$pseq->display_id('hello');

if( $pseq->seq ne 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' ) {
	print "not ok 2\n";
	print STDERR "Sequence is ",$pseq->seq,"\n";

} else {
	print "ok 2\n";
}

if( $pseq->subseq(3,7) ne 'GGGGT' ) {
	print "not ok 3\n";
	print STDERR "Subseq is ",$pseq->subseq(3,7),"\n";
} else {
	print "ok 3\n";
}











