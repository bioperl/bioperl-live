# -*-Perl-*-
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
BEGIN { $| = 1; print "1..7\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}
$loaded = 1;
print "ok 1\n";
use lib '../';

use Bio::SeqIO;

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my $seqio = new Bio::SeqIO(-format=>'largefasta',
			   -file=>'t/genomic-seq.fasta');

test 2, defined $seqio, 'cannot instantiate Bio::SeqIO::largefasta';

my $pseq = $seqio->next_seq();
$pseq->moltype('dna');
$pseq->desc('this is my description');


test 3, defined $pseq, 'could not call next_seq';
test 4, $pseq->length() > 0, "could not call length, seq was empty";
test 5, length($pseq->subseq(100, 299)) == 200, 'error in subseq'; 
test 6, $pseq->trunc(100,199)->length() == 100, 'error in trunc'; 
test 7, $pseq->moltype() eq 'dna', 'moltype was ' . $pseq->moltype();
test 8, $pseq->display_id(), "no display id";
test 9, $pseq->accession_number(), "no accession";
test 10, $pseq->desc, "no description ";
