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
    plan tests => 34;
}

use Bio::Tools::Genscan;
use Bio::Tools::Genemark;
use Bio::Tools::MZEF;  
use Bio::SeqIO;
use Bio::Root::IO;

# Genscan report
my $genscan = Bio::Tools::Genscan->new('-file' => Bio::Root::IO->catfile("t","data","genomic-seq.genscan"));
ok $genscan;

# original sequence
my $seqin = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","genomic-seq.fasta"),
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
	$fea = ($gene->exons_ordered())[0];
	ok $fea->abs_strand(), -1, 
	     "strand mismatch (".$fea->abs_strand()." instead of -1)";
	$fea = ($gene->poly_A_site());
	ok $fea->score(), 1.05, 
             "score mismatch (".$fea->score()." instead of 1.05)";
    }
    if($pred_num == 2) {
	$fea = ($gene->exons_ordered("Initial"))[0];
	ok $fea->abs_strand(), 1, 
	"strand mismatch (".$fea->abs_strand()." instead of 1)";
	ok $fea->score(), 4.46, 
             "score mismatch (".$fea->score()." instead of 4.46)";
    }
    if($pred_num == 3) {
	my @exons = $gene->exons("Initial");
	ok scalar(@exons), 0, 
	     "initial exons (".scalar(@exons)." instead of 0)";
	$fea = ($gene->exons_ordered())[-1];
	ok $fea->score(),  1.74, 
             "score mismatch (".$fea->score()." instead of 1.74)";
    }
    if($seq) {
	$prtseq = $gene->predicted_protein()->seq();
        $cds = $gene->cds();
	ok($cds) || print STDERR "# no CDS for prediction $pred_num; protein: $prtseq\n";
	$tr_cds = $cds->translate()->seq();
	$tr_cds =~ s/\*$//;
	ok( lc($prtseq), lc($tr_cds),
	    "predicted and extracted protein seqs don't match for prediction $pred_num, gene $gene (RelRangeString is ".$gene->toRelRangeString( 'both', 'plus' ).".");
    }
}

# Genscan report with no genes predicted
my $null_genscan = Bio::Tools::Genscan->new('-file' => Bio::Root::IO->catfile("t","data","no-genes.genscan"));
ok $null_genscan;
my $no_gene = $null_genscan->next_prediction;
my @exons = $no_gene->exons;
ok($#exons,-1);

# MZEF report
my $mzef = Bio::Tools::MZEF->new('-file' => Bio::Root::IO->catfile("t","data","genomic-seq.mzef"));
ok $mzef;

my $exon_num = 0;
my $gene = $mzef->next_prediction();

ok($gene->exons, 23);

# Genemark testing:
my $genemark = Bio::Tools::Genemark->new('-file' => Bio::Root::IO->catfile("t", "data", "genemark.out"));

my $gmgene = $genemark->next_prediction();
ok $gmgene->seq_id(), "Hvrn.contig8";
ok $genemark->analysis_date(), "Thu Mar 22 10:25:00 2001";

my $i = 0;
my @num_exons = (1,5,2,1,9,5,3,2,3,2,1,2,7);
while($gmgene = $genemark->next_prediction()) {
    $i++;
    my @gmexons = $gmgene->exons_ordered();
    ok scalar(@gmexons), $num_exons[$i];

    if($i == 5) {
	my $gmstart = $gmexons[-1]->abs_start( 'plus' );
	ok $gmstart, 23000;

	my $gmend = $gmexons[-1]->abs_end( 'plus' );
	ok $gmend, 23061;
    }
}

