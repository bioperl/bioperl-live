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
BEGIN { $| = 1; print "1..24\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::BPlite;

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

open FH, "t/blast.report";
my $report = Bio::Tools::BPlite->new(-fh=>\*FH);
test 2, $report;
my $sbjct = $report->nextSbjct;
test 3, $sbjct;
my $hsp = $sbjct->nextHSP;
test 4, $hsp;

test 5, $report->query eq "gi|1401126 (504 letters) ";
test 6, $report->database eq 'Non-redundant GenBank+EMBL+DDBJ+PDB sequences';
test 7, $sbjct->name eq 'gb|U49928|HSU49928 Homo sapiens TAK1 binding protein (TAB1) mRNA, complete cds. ';
test 8, $hsp->bits == 1009;
test 9, $hsp->score == 2580;
test 10, $hsp->percent == 100;
test 11, $hsp->P == 0.0;
test 12, $hsp->match == 504;
test 13, $hsp->positive == 504;
test 14, $hsp->length == 504;
test 15, $hsp->querySeq eq $seq;
test 16, $hsp->sbjctSeq eq $seq;
test 17, $hsp->homologySeq eq $seq;
test 18, $hsp->query->start == 1;
test 19, $hsp->query->end == 504;
test 20, $hsp->query->seqname eq $report->query;
test 21, $hsp->query->primary_tag eq "similarity";
test 22, $hsp->query->source_tag eq "BLAST";
test 23, $hsp->subject->length == 1512;

close FH;

test 24, "everything fine";





