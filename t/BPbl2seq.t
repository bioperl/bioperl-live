# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# $Id$

use Test;
use strict;
BEGIN { plan tests => 18 }

use Bio::Tools::BPbl2seq;
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

open FH, "t/bl2seq.out";
my $report = Bio::Tools::BPbl2seq->new(\*FH);
ok ref($report), qr/Bio::Tools::BPbl2seq/, " no report";
ok defined($report->query),1, " no query";
ok $report->score, 481, "wrong score";
ok $report->bits, 191, "wrong score in bits ";
ok $report->percent, 35.1, "wrong match percent";
ok $report->P, 2e-53, "wrong expectation value ";
ok $report->match, 111, "wrong number of matches ";
ok $report->positive, 167, "wrong number of positives";
ok $report->length, 316, "wrong length";
ok $report->querySeq, qr/QFL/ , "bad query sequence";
ok $report->sbjctSeq, qr/RFAR/ , "bad hit sequence";
ok $report->homologySeq, qr/PVKN/ , "bad homology sequence";
ok $report->query->start, 28, "wrong query start";
ok $report->query->end, 343, "wrong query end";
ok $report->subject->start, 60, "wrong hit start ";
ok $report->subject->end, 360, "wrong hit end";
ok $report->subject->seqname, qr/ALEU_HORVU/ , "wrong hit name";
close FH;
