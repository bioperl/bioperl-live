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
BEGIN { $| = 1; print "1..13\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::Genscan;
use Bio::Tools::MZEF;  # THIS IS STILL TODO!
use Bio::SeqIO;

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

# Genscan report
my $genscan = Bio::Tools::Genscan->new('-file' => "t/genomic-seq.genscan");
test 2, $genscan, "unable to open test output file";

# original sequence
my $seqin = Bio::SeqIO->new('-file' => "t/genomic-seq.fasta",
			    '-format' => "fasta");
test 3, $seqin, "unable to open test sequence file";
$seq = $seqin->next_seq();
$seqin->close();
test 4, $seq, "unable to read seq from seq file";

# scan through the report
my $fea;
my $pred_num = 0;
my ($prtseq, $cds, $tr_cds);
while($gene = $genscan->next_prediction()) {
    $gene->attach_seq($seq) if $seq;
    $pred_num++;

    if($pred_num == 1) {
	$fea = ($gene->exons())[0];
	test 4+(($pred_num-1)*3)+1, $fea->strand() == -1, 
	     "strand mismatch (".$fea->strand()." instead of -1)";
	$fea = ($gene->poly_A_sites())[0];
	test 4+(($pred_num-1)*3)+2, $fea->score() == 1.05, 
             "score mismatch (".$fea->score()." instead of 1.05)";
    }
    if($pred_num == 2) {
	$fea = ($gene->exons("Initial"))[0];
	test 4+(($pred_num-1)*3)+1, $fea->strand() == 1, 
	     "strand mismatch (".$fea->strand()." instead of 1)";
	test 4+(($pred_num-1)*3)+2, $fea->score() == 4.46, 
             "score mismatch (".$fea->score()." instead of 4.46)";
    }
    if($pred_num == 3) {
	my @exons = $gene->exons("Initial");
	test 4+(($pred_num-1)*3)+1, scalar(@exons) == 0, 
	     "initial exons (".scalar(@exons)." instead of 0)";
	$fea = ($gene->exons())[0];
	test 4+(($pred_num-1)*3)+2, $fea->score() == 1.74, 
             "score mismatch (".$fea->score()." instead of 1.74)";
    }
    if($seq) {
	$prtseq = $gene->predicted_protein()->seq();
        $cds = Bio::PrimarySeq->new('-seq' => $gene->cds(1),
                                    '-id' => "cds_" . $gene->primary_tag());
	$tr_cds = $cds->translate()->seq();
	$tr_cds =~ s/\*$//;
	test 4+(($pred_num-1)*3)+3, lc($prtseq) eq lc($tr_cds), "predicted and extracted protein seqs don't match";
    }
}


