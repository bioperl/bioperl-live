# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
    plan tests => 76;
}

use Bio::Tools::BPlite;
use Bio::Root::IO;
ok(1);

my $seq =
    "MAAQRRSLLQSEQQPSWTDDLPLCHLSGVGSASNRSYSADGKGTESHPPEDSWLKFRSENN".
    "CFLYGVFNGYDGNRVTNFVAQRLSAELLLGQLNAEHAEADVRRVLLQAFDVVERSFLESID".
    "DALAEKASLQSQLPEGVPQHQLPPQYQKILERLKTLEREISGGAMAVVAVLLNNKLYVANV".
    "GTNRALLCKSTVDGLQVTQLNVDHTTENEDELFRLSQLGLDAGKIKQVGIICGQESTRRIG".
    "DYKVKYGYTDIDLLSAAKSKPIIAEPEIHGAQPLDGVTGFLVLMSEGLYKALEAAHGPGQA".
    "NQEIAAMIDTEFAKQTSLDAVAQAVVDRVKRIHSDTFASGGERARFCPRHEDMTLLVRNFG".
    "YPLGEMSQPTPSPAPAAGGRVYPVSVPYSSAQSTSKTSVTLSLVMPSQGQMVNGAHSASTL".
    "DEATPTLTNQSPTLTLQSTNTHTQSSSSSSDGGLFRSRPAHSLPPGEDGRVEPYVDFAEFY".
    "RLWSVDHGEQSVVTAP";

open FH, Bio::Root::IO->catfile("t","data","blast.report");
my $report = Bio::Tools::BPlite->new(-fh=>\*FH);
ok $report->isa('Bio::Tools::BPlite');
my $sbjct = $report->nextSbjct;
ok defined $sbjct;
my $hsp = $sbjct->nextHSP;
ok defined $hsp;

ok $report->query, "gi|1401126 (504 letters) ";
ok $report->database, 'Non-redundant GenBank+EMBL+DDBJ+PDB sequences';
ok $sbjct->name, 'gb|U49928|HSU49928 Homo sapiens TAK1 binding protein (TAB1) mRNA, complete cds. ';
ok $hsp->bits, 1009;
ok $hsp->score, 2580;
ok $hsp->percent, 100;
ok $hsp->P, '0.0';
ok $hsp->match, 504;
ok $hsp->positive, 504;
ok $hsp->length, 504;
ok $hsp->querySeq, $seq;
ok $hsp->sbjctSeq, $seq;
ok $hsp->homologySeq, $seq;
ok $hsp->query->start, 1;
ok $hsp->query->end, 504;
ok $hsp->query->seqname, $report->query;
ok $hsp->query->primary_tag, "similarity";
ok $hsp->query->source_tag, "BLAST";
ok $hsp->subject->length, 1512;

close FH;

# Verify that BPlite is properly parsing PHIBLAST reports as well

my $report2 = Bio::Tools::BPlite->new(-file=>Bio::Root::IO->catfile("t","data","phi.out"));

ok $report2->pattern, "P-E-E-Q";
ok $report2->query_pattern_location->[0], 23;
ok $report2->query_pattern_location->[1], 120;
my $sbjct2 = $report2->nextSbjct;
ok $sbjct2->name =~ /4988/;
my $hsp2 = $sbjct2->nextHSP;
ok $hsp2->subject->end, 343;

close FH;


# test SeqAnalysisParserI

# tests 29-38 are just counting to see that we get the expected number 
# of features
my $parser = new Bio::Tools::BPlite(-file => Bio::Root::IO->catfile("t","data","blast.report"));
while( $parser->next_feature ) {
    ok(1);
}

$parser = new Bio::Tools::BPlite(-file => Bio::Root::IO->catfile("t","data","cysprot.tblastn"));
while( $parser->next_feature ) {
    ok(1);
}

$parser = new Bio::Tools::BPlite(-file => Bio::Root::IO->catfile("t",
								 "data",
								 "short.blx"));
ok($parser);
my $count = 0;
while( $parser->next_feature ) {
    ok(1);    
}


