## Bioperl Test Harness Script for Modules
##


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


use Bio::SimpleAlign;
use Bio::AlignIO;

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


## Now we test Bio::AlignIO::stockholm input
$str = Bio::AlignIO->new(-file=> 't/testaln.stockholm','-format' => 'stockholm');
$aln = $str->next_aln();
test 2, $aln->{order}->{'0'} eq '1433_LYCES-9-246', " failed stockholm format test";

## Now we test Bio::AlignIO::pfam

$str = Bio::AlignIO->new(-file=> 't/testaln.pfam');
$aln = $str->next_aln();
test 3,$aln->{order}->{'0'} eq '1433_LYCES-9-246', " failed pfam input test";

$strout = Bio::AlignIO->new(-file=> '>t/testout.pfam', '-format' => 'pfam');
$status = $strout->write_aln($aln);
test 4,$status, " failed pfam output test";


## Now we test Bio::AlignIO::msf

$str = Bio::AlignIO->new(-file=> 't/testaln.msf');
$aln = $str->next_aln();
test 5, $aln->{order}->{'0'} eq  '1433_LYCES-9-246', " failed msf input test";



$strout = Bio::AlignIO->new(-file=> '>t/testout.msf', '-format' => 'msf');
$status = $strout->write_aln($aln);
test 6, $status, "  failed msf output test";


## Now we test Bio::AlignIO::fasta

$str = Bio::AlignIO->new(-file=> 't/testaln.fasta', '-format' => 'fasta');
$aln = $str->next_aln();
test 7, $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378', " failed fasta input test ";


$strout = Bio::AlignIO->new(-file=> '>t/testout.fasta', '-format' => 'fasta');
$status = $strout->write_aln($aln);
test 8, $status, "  failed fasta output test";

## Now we test Bio::AlignIO::selex

$str = Bio::AlignIO->new(-file=> 't/testaln.selex','-format' => 'selex');
$aln = $str->next_aln();
test 9, $aln->{order}->{'0'} eq 'AK1H_ECOLI-114-431', " failed selex format test ";


$strout = Bio::AlignIO->new(-file=> '>t/testout.selex', '-format' => 'selex');
$status = $strout->write_aln($aln);
test 10, $status, "  failed selex output test";

## Now we test Bio::AlignIO::mase input
$str = Bio::AlignIO->new(-file=> 't/testaln.mase','-format' => 'mase');
$aln = $str->next_aln();
test 11, $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378', " failed mase input test ";

## Now we test Bio::AlignIO::prodom input
$str = Bio::AlignIO->new(-file=> 't/testaln.prodom','-format' => 'prodom');
$aln = $str->next_aln();
test 12, $aln->{order}->{'0'} eq 'P04777-1-33', " failed prodom input test ";

## Now we test Bio::AlignIO::clustalw output writing
$strout = Bio::AlignIO->new(-file=> '>t/testaln.clustal', '-format' => 'clustalw');
$status = $strout->write_aln($aln);
test 13, $status, "  failed clustalw (.aln) output test";


# Testing filehandle manipulations

my $in  = Bio::AlignIO->newFh(-file => "t/testaln.fasta", '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(-file => ">t/testout2.pfam", '-format' => 'pfam');
while ( $aln = <$in>) {
	test 14, $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378', "  failed filehandle input test  ";
	$status = print $out $aln;
}
test 15, $status, "  failed filehandle output test";

