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
BEGIN { $| = 1; print "1..16\n"; 
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

my $tmpfile = 't/largefastatest.out';
my $seqio = new Bio::SeqIO('-format'=>'largefasta',
			   '-file'  =>'t/genomic-seq.fasta');
test 2, defined $seqio, 'cannot instantiate Bio::SeqIO::largefasta';

my $pseq = $seqio->next_seq();
$pseq->moltype('dna');
$pseq->desc('this is my description');
my $plength = $pseq->length();
my $last_3 = $pseq->subseq($plength-3,$plength);

test 3, defined $pseq, 'could not call next_seq';
test 4, $plength > 0, "could not call length, seq was empty";
test 5, length($pseq->subseq(100, 299)) == 200, 'error in subseq'; 
test 6, $pseq->trunc(100,199)->length() == 100, 'error in trunc'; 
test 7, $pseq->moltype() eq 'dna', 'moltype was ' . $pseq->moltype();
test 8, $pseq->display_id(), "no display id";
test 9, $pseq->accession_number(), "no accession";
test 10, $pseq->desc, "no description ";

test 11, open(OUT, ">$tmpfile"), 'could not open output file';

my $seqout = new Bio::SeqIO('-format' => 'largefasta',
			    '-fh'     => \*OUT );
test 12, defined $seqout, 'could not open seq with outputstream';

test 13, $seqout->write_seq($pseq), 'could not write seq';
$seqout->close();
close(OUT);
$seqin = new Bio::SeqIO('-format' => 'largefasta',
			'-file'   => $tmpfile);
my $pseq2 = $seqin->next_seq;
test 14, ( $plength == $pseq2->length() ), "written file was not same length as expected";
test 15, ( $pseq->display_id() eq $pseq2->display_id() ), "display ids were not identical as expected";
test 16, ( $pseq->desc() eq $pseq2->desc() ), "description was not identical (" . $pseq->desc() . "," . $pseq2->desc() . ")";

END {
    unlink $tmpfile;
}
