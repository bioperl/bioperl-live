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
    plan tests => 61; 
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
ok int $hsp->percent, 34, "wrong match percent";
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
ok $hsp->hit->start, 60, "wrong hit start ";
ok $hsp->hit->end, 360, "wrong hit end";
ok $report->sbjctName =~ /ALEU_HORVU/;# "wrong hit name";

$report = new Bio::Tools::BPbl2seq(-file => Bio::Root::IO->catfile("t","data","bl2seq.bug940.out"));

$report->verbose(2);
ok $report->isa('Bio::Tools::BPbl2seq');# " no report";
ok defined($report->sbjctName),1, " no query";
$hsp = $report->next_feature;
ok $hsp->score, 1626, "wrong score";
ok $hsp->bits, 637, "wrong score in bits ";
ok int $hsp->percent, 66, "wrong match percent";
ok $hsp->P == 0.0;# "wrong expectation value ";
ok $hsp->match, 311, "wrong number of matches ";
ok $hsp->positive, 330, "wrong number of positives";
ok $hsp->start, 121, 'wrong starting position';
ok $hsp->end, 469, 'wrong ending position';
ok $hsp->length, 349, "wrong length";
ok $hsp->querySeq =~ /^MGN/; #"bad query sequence";
ok $hsp->sbjctSeq =~ /^MGN/;#"bad hit sequence";
ok $hsp->homologySeq =~ /^MGN/;# , "bad homology sequence";
ok $hsp->query->start, 121, "wrong query start";
ok $hsp->query->end, 469, "wrong query end";
ok $hsp->hit->start, 1, "wrong hit start ";
ok $hsp->hit->end, 469, "wrong hit end";
ok $hsp->hit->seqname =~ /gi|4507985/;# "wrong hit name";
ok $hsp->gaps, 120;

$hsp = $report->next_feature;
ok($hsp);
ok $hsp->score, 1524, "wrong score";
ok $hsp->bits, 598, "wrong score in bits ";
ok int $hsp->percent, 61, "wrong match percent";
ok $hsp->P == 1e-175, 1,"wrong expectation value ";
ok $hsp->match, 275, "wrong number of matches ";
ok $hsp->positive, 324, "wrong number of positives";
ok $hsp->start, 6, 'wrong starting position';
ok $hsp->end, 420, 'wrong ending position';
ok $hsp->length, 415, "wrong length";
ok $hsp->querySeq =~ /^EKPY/; #"bad query sequence";
ok $hsp->sbjctSeq =~ /^EKPY/;#"bad hit sequence";
ok $hsp->homologySeq =~ /^EKPY/;# , "bad homology sequence";
ok $hsp->query->start, 6, "wrong query start";
ok $hsp->query->end, 420, "wrong query end";
ok $hsp->hit->start, 22, "wrong hit start ";
ok $hsp->hit->end, 464, "wrong hit end";
ok $hsp->hit->seqname =~ /gi|4507985/;# "wrong hit name";
ok $hsp->gaps, 30;

$report = new Bio::Tools::BPbl2seq(-file =>  Bio::Root::IO->catfile("t","data","empty.bl2seq"));
$report->verbose(2);
ok( $report->isa('Bio::Tools::BPbl2seq'),1, " no report found");
ok $report->sbjctName, '',"query found where none expected";
