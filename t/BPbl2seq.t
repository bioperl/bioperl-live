## Bioperl Test Harness Script for Modules
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
BEGIN { $| = 1; print "1..18\n";
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::BPbl2seq;

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
open FH, "t/bl2seq.out";
my $report = Bio::Tools::BPbl2seq->new(\*FH);
test 2, $report, " no report";
test 3, $report->query, " no query";
test 4, $report->score == 481, "wrong score";
test 5, $report->bits == 191, "wrong score in bits ";
test 6, $report->percent == 35.1, "wrong match percent";
test 7, $report->P == 2e-53, "wrong expectation value ";
test 8, $report->match == 111, "wrong number of matches ";
test 9, $report->positive == 167, "wrong number of positives";
test 10, $report->length == 316, "wrong length";
test 11, $report->querySeq =~ /QFL/ , "bad query sequence";
test 12, $report->sbjctSeq =~ /RFAR/ , "bad hit sequence";
test 13,  $report->homologySeq =~ /PVKN/ , "bad homology sequence";
test 14, $report->query->start == 28, "wrong query start";
test 15, $report->query->end == 343, "wrong query end";
test 16, $report->subject->start == 60, "wrong hit start ";
test 17, $report->subject->end == 360, "wrong hit end";
test 18, $report->subject->seqname =~ /ALEU_HORVU/ , "wrong hit name";
close FH;
__END__
