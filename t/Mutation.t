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
BEGIN { $| = 1; print "1..15\n";
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::LiveSeq::Mutation;

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

my $a = Bio::LiveSeq::Mutation->new();
test 2, $a;
test 3, $a->seq('aaa');
test 4, $a->seqori('ggg');
test 5, $a->pos(-4);
test 6, $a->pos(5);
test 7, ($a->len == 3);
test 8, $a->len(9);
test 9, $a->transpos(55);
test 10, $a->issue(1);
test 11, $a->label(57);
test 12, $a->prelabel(57);
test 13, $a->postlabel(57);
test 14, $a->lastlabel(57);

#constuctor test
$b = Bio::LiveSeq::Mutation->new(-seq=>'AC',
				 -seqori => 'GG'
				 -pos => 5,
				 -len => 2
				 );
test 15,  $b;



