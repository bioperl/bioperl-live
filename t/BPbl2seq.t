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
#use Bio::Tools::BPlite;
use Bio::Root::IO;
ok(1);

my $report = new Bio::Tools::BPbl2seq(-file => Bio::Root::IO->catfile("t","data","bl2seq.out"));
$report->verbose(2);
ok $report->isa('Bio::Tools::BPbl2seq');# " no report";
ok defined($report->sbjctName),1, " no hit";
my $hsp = $report->next_feature;
ok $hsp->score, 481, "wrong score";
ok $hsp->bits, 191, "wrong score in bits ";
ok $hsp->percent, 35.1, "wrong match percent";
ok $hsp->P == 2e-53;# "wrong expectation value ";
ok $hsp->match, 111, "wrong number of matches ";
ok $hsp->positive, 167, "wrong number of positives";
ok $hsp->start, 28, 'wrong starting position';
ok $hsp->end, 343, 'wrong ending position';
ok $hsp->length, 316, "wrong length";
ok $hsp->querySeq =~ /QFL/; #"bad query sequence";
ok $hsp->sbjctSeq =~ /RFAR/;#"bad hit sequence";
ok $hsp->homologySeq =~ /PVKN/;# , "bad homology sequence";
ok $hsp->query->start, 28, "wrong query start";
ok $hsp->query->end, 343, "wrong query end";
ok $hsp->subject->start, 60, "wrong hit start ";
ok $hsp->subject->end, 360, "wrong hit end";
ok $report->sbjctName =~ /ALEU_HORVU/;# "wrong hit name";
