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
    plan tests => 14;
}

use Bio::Tools::Genscan;
use Bio::Tools::MZEF;  # THIS IS STILL TODO!
use Bio::SeqIO;

# Genscan report
my $genscan = Bio::Tools::Genscan->new('-file' => "t/genomic-seq.genscan");
ok $genscan;

# original sequence
my $seqin = Bio::SeqIO->new('-file' => "t/genomic-seq.fasta",
			    '-format' => "fasta");
ok $seqin;
my $seq = $seqin->next_seq();
$seqin->close();
ok $seq;

# scan through the report
my $fea;
my $pred_num = 0;
my ($prtseq, $cds, $tr_cds);
while(my $gene = $genscan->next_prediction()) {
    $gene->attach_seq($seq) if $seq;
    $pred_num++;

    if($pred_num == 1) {
	$fea = ($gene->exons())[0];
	ok $fea->strand(), -1, 
	     "strand mismatch (".$fea->strand()." instead of -1)";
	$fea = ($gene->poly_A_sites())[0];
	ok $fea->score(), 1.05, 
             "score mismatch (".$fea->score()." instead of 1.05)";
    }
    if($pred_num == 2) {
	$fea = ($gene->exons("Initial"))[0];
	ok $fea->strand(), 1, 
	"strand mismatch (".$fea->strand()." instead of 1)";
	ok $fea->score(), 4.46, 
             "score mismatch (".$fea->score()." instead of 4.46)";
    }
    if($pred_num == 3) {
	my @exons = $gene->exons("Initial");
	ok scalar(@exons), 0, 
	     "initial exons (".scalar(@exons)." instead of 0)";
	$fea = ($gene->exons())[0];
	ok $fea->score(),  1.74, 
             "score mismatch (".$fea->score()." instead of 1.74)";
    }
    if($seq) {
	$prtseq = $gene->predicted_protein()->seq();
        $cds = Bio::PrimarySeq->new('-seq' => $gene->cds(1),
                                    '-id' => "cds_" . $gene->primary_tag());
	$tr_cds = $cds->translate()->seq();
	$tr_cds =~ s/\*$//;
	ok( lc($prtseq), lc($tr_cds),
	    "predicted and extracted protein seqs don't match");
    }
}

# MZEF report
my $mzef = Bio::Tools::MZEF->new('-file' => "t/genomic-seq.mzef");
ok $mzef;

my $exon_num = 0;
my $gene = $mzef->next_prediction();

ok($gene->exons, 23);
