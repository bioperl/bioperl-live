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
BEGIN { $| = 1; print "1..11\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Seq;
use Bio::SeqFeature::Generic;

$loaded = 1;
print "ok 1\n";    # 1st test passes.


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    $msg = '' if ( !defined $msg );
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACT',
                        -desc=>'Sample Bio::Seq object',
			-moltype => 'dna' );
test 2, 1;

$trunc = $seq->trunc(1,4);

test 3, 1;

test 4, ( $trunc->seq() eq 'ACTG' );

$trans = $seq->translate();

test 5, ( $trans->seq() eq 'TVAST' );

# test ability to get str function

$t = $seq->seq();
test 6, ( $t eq 'ACTGTGGCGTCAACT' );

$seq = Bio::Seq->new(-seq=>'actgtggcgtcaact',
		     -desc=>'Sample Bio::Seq object',
		     -display_id => 'something',
		     -accession_number => 'accnum',
		     -moltype => 'dna' );
test 7, defined $seq;


$trans = $seq->translate();

test 8, ( $trans->seq() eq 'TVAST' );

# basic methods

test 9, ( $seq->id() eq 'something' && 
	$seq->accession_number eq 'accnum' ), "saw ".$seq->id.":".
	$seq->accession_number.":".$seq->primary_id;

my $subseq = $seq->subseq(5, 9);

test 10, ( $seq->subseq(5, 9) eq 'tggcg'), "subseq(5,9) was ".
    $seq->subseq(5,9). " when I expected tggcg\n";

my $newfeat = Bio::SeqFeature::Generic->new( -start => 10,
					     -end => 12,
						-primary => 'silly',
						-source => 'stuff');


$seq->add_SeqFeature($newfeat);
test 11, 1;
