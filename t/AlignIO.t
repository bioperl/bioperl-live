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


## Now we test Bio::AlignIO::stockholm input
$str = Bio::AlignIO->new(-file=> 't/testaln.stockholm','-format' => 'stockholm');
$aln = $str->next_aln();

if( $aln->{order}->{'0'} eq '1433_LYCES-9-246' ) {
    print "ok 2\n";
} else {
    print "not ok 2 , failed stockholm format test \n";	
}


## Now we test Bio::AlignIO::pfam

$str = Bio::AlignIO->new(-file=> 't/testaln.pfam');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq '1433_LYCES-9-246' ) {
    print "ok 3\n";
} else {
    print "not ok 3, failed pfam input test\n";	
}
$strout = Bio::AlignIO->new(-file=> '>t/testout.pfam', '-format' => 'pfam');
$status = $strout->write_aln($aln);
if($status) {
    print "ok 4\n";
} else {
    print "not ok 4, failed pfam output test\n";	
}

## Now we test Bio::AlignIO::msf

$str = Bio::AlignIO->new(-file=> 't/testaln.msf');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq  '1433_LYCES-9-246' ) {
    print "ok 5\n";
} else {
    print "not ok 5, failed msf input test\n";	
}




$strout = Bio::AlignIO->new(-file=> '>t/testout.msf', '-format' => 'msf');
$status = $strout->write_aln($aln);
if($status) {
    print "ok 6\n";
} else {
    print "not ok 6, failed msf output test\n";	
}

## Now we test Bio::AlignIO::fasta

$str = Bio::AlignIO->new(-file=> 't/testaln.fasta', '-format' => 'fasta');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378' ) {
    print "ok 7\n";
} else {
    print "not ok 7, failed fasta input test \5n";	
}

$strout = Bio::AlignIO->new(-file=> '>t/testout.fasta', '-format' => 'fasta');
$status = $strout->write_aln($aln);
if($status) {
    print "ok 8\n";
} else {
    print "not ok 8, failed fasta output test\n";	
}

## Now we test Bio::AlignIO::selex

$str = Bio::AlignIO->new(-file=> 't/testaln.selex','-format' => 'selex');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq 'AK1H_ECOLI-114-431' ) {
    print "ok 9\n";
} else {
    print "not ok 9, failed selex format test\n";	
}

$strout = Bio::AlignIO->new(-file=> '>t/testout.selex', '-format' => 'selex');
$status = $strout->write_aln($aln);
if($status) {
    print "ok 10\n";
} else {
    print "not ok 10, failed selex output test\n";	
}


## Now we test Bio::AlignIO::mase input
$str = Bio::AlignIO->new(-file=> 't/testaln.mase','-format' => 'mase');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378' ) {
    print "ok 11\n";
} else {
    print "not ok 11, failed mase input test\n";	
}


## Now we test Bio::AlignIO::prodom input
$str = Bio::AlignIO->new(-file=> 't/testaln.prodom','-format' => 'prodom');
$aln = $str->next_aln();
if( $aln->{order}->{'0'} eq 'P04777-1-33' ) {
    print "ok 12\n";
} else {
    print "not ok 12, failed prodom input test\n";	
}


## Now we test Bio::AlignIO::clustalw output writing
$strout = Bio::AlignIO->new(-file=> '>t/testaln.clustal', '-format' => 'clustalw');
$status = $strout->write_aln($aln);
if($status) {
    print "ok 13\n";
} else {
    print "not ok 13, failed clustalw (.aln) output test\n";	
}



# Testing filehandle manipulations

    my $in  = Bio::AlignIO->newFh(-file => "t/testaln.fasta", '-format' => 'fasta');
    my $out = Bio::AlignIO->newFh(-file => ">t/testout2.pfam", '-format' => 'pfam');
     while ( $aln = <$in>) {
	if( $aln->{order}->{'0'} eq 'AK1H_ECOLI-1-378' ) {
    		print "ok 14\n";
	} else {
    		print "not ok 14, failed filehandle input test \n";	
	}
	$status = print $out $aln;
	if($status) {
	    print "ok 15\n";
	} else {
	    print "not ok 15, failed filehandle output test\n";	
	}
    }
