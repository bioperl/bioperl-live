# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# $Id$

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 20; 
}

use Bio::Tools::BPbl2seq;
use Bio::Root::IO;
ok(1);

my $report = new Bio::Tools::BPbl2seq(-file => Bio::Root::IO->catfile("t","bl2seq.out"));
$report->verbose(2);
ok $report->isa('Bio::Tools::BPbl2seq');# " no report";
ok defined($report->query),1, " no query";
ok $report->score, 481, "wrong score";
ok $report->bits, 191, "wrong score in bits ";
ok $report->percent, 35.1, "wrong match percent";
ok $report->P , 2e-53;# "wrong expectation value ";
ok $report->match, 111, "wrong number of matches ";
ok $report->positive, 167, "wrong number of positives";
ok $report->start, 28, 'wrong starting position';
ok $report->end, 343, 'wrong ending position';
ok $report->length, 316, "wrong length";
ok $report->querySeq =~ /QFL/; #"bad query sequence";
ok $report->sbjctSeq =~ /RFAR/;#"bad hit sequence";
ok $report->homologySeq =~ /PVKN/;# , "bad homology sequence";
ok $report->query->start, 28, "wrong query start";
ok $report->query->end, 343, "wrong query end";
ok $report->subject->start, 60, "wrong hit start ";
ok $report->subject->end, 360, "wrong hit end";
ok $report->subject->seqname =~ /ALEU_HORVU/;# "wrong hit name";
