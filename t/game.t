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
BEGIN { $| = 1; print "1..8\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$str = Bio::SeqIO->new(-file=> 't/test.game', '-format' => 'game');
test 2, ($str);

test 3, ($seq = $str->next_primary_seq()), 'failed to read game primary_seq from stream'; 

test 4, ( $seq->id eq 'AE003417' );

$str2 = Bio::SeqIO->new(-file=> 't/test.game', '-format' => 'game');
test 5, ($str2);

test 6, ($seq = $str2->next_seq()), 'failed to read game seq from stream'; 

#$str2->write_seq($seq);

test 7, ( $seq->id eq 'AE003417' );

@feats = $seq->all_SeqFeatures();

test 8, ( @feats[0]->primary_tag eq 'exon' );

