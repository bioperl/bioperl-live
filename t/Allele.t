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
BEGIN { $| = 1; print "1..13\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Variation::Allele;

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

my $a = Bio::Variation::Allele->new(-seq=>'ACTGACTGACTG',
			-display_id => 'new-id',
			-moltype => 'dna',
			-accession_number => 'X677667',
                        -desc=>'Sample Bio::Seq object');
test 2, defined $a && ref($a) =~ /Bio::Variation::Allele/;


test 3, $a->accession_number() && $a->seq() && $a->display_id();


$trunc = $a->trunc(1,4);

test 4, defined $trunc;

test 5, ( $trunc->seq() eq 'ACTG' ), "Expecting ACTG. Got ". $trunc->seq();

$rev = $a->revcom();

test 6, defined $rev;

test 7, ( $rev->seq() eq 'CAGTCAGTCAGT' );

$a->is_reference(1);

test 8, defined $a;

test 9, ( $a->is_reference );

$a->repeat_unit('ACTG');
test 10, defined $a;

test 11, ( $a->repeat_unit eq 'ACTG' );

$a->repeat_count(3);
test 12, defined $a;

test 13, ( $a->repeat_count == 3 );

