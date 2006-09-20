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
    plan tests => 79;
}

use Bio::Tools::Genscan;
use Bio::Tools::Genemark;
use Bio::Tools::Glimmer;
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
	$fea = ($gene->exons())[0];
	ok $fea->strand(), -1, 
	     "strand mismatch (".$fea->strand()." instead of -1)";
	$fea = ($gene->poly_A_site());
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
        $cds = $gene->cds();
	ok($cds) || print STDERR "# no CDS for prediction $pred_num; protein: $prtseq\n";
	$tr_cds = $cds->translate()->seq();
	$tr_cds =~ s/\*$//;
	ok( lc($prtseq), lc($tr_cds),
	    "predicted and extracted protein seqs don't match");
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
my $genemark = Bio::Tools::Genemark->new('-file' => Bio::Root::IO->catfile(qw(t data genemark.out)));

my $gmgene = $genemark->next_prediction();
ok $gmgene->seq_id(), "Hvrn.contig8";
ok $genemark->analysis_date(), "Thu Mar 22 10:25:00 2001";

my $i = 0;
my @num_exons = (1,5,2,1,9,5,3,2,3,2,1,2,7);
while($gmgene = $genemark->next_prediction()) {
    $i++;
    my @gmexons = $gmgene->exons();
    ok scalar(@gmexons), $num_exons[$i];

    if($i == 5) {
	my $gmstart = $gmexons[0]->start();
	ok $gmstart, 23000;

	my $gmend = $gmexons[0]->end();
	ok $gmend, 23061;
    }
}

# Glimmer testing (GlimmerM)
my $glimmer = new Bio::Tools::Glimmer('-file' => Bio::Root::IO->catfile(qw(t data glimmer.out)));
my $glimmergene = $glimmer->next_prediction;

ok($glimmergene);
ok($glimmergene->seq_id, 'BAC1Contig11');
ok($glimmergene->source_tag, 'GlimmerM_3.0');
ok($glimmergene->primary_tag, 'transcript');
ok(($glimmergene->get_tag_values('Group'))[0], 'GenePrediction1');
my @glim_exons = $glimmergene->exons;
ok(scalar (@glim_exons), 5);
ok($glim_exons[0]->start, 13907);
ok($glim_exons[0]->end, 13985);
ok($glim_exons[0]->strand, 1);
ok(($glim_exons[0]->get_tag_values('Group'))[0], 'GenePrediction1');

@num_exons = (0,5,3, 1, 6, 3);
$i = 1;
while($glimmergene = $glimmer->next_prediction()) {
    $i++;
    ok(($glimmergene->get_tag_values('Group'))[0],"GenePrediction$i");
    @glim_exons = $glimmergene->exons();    
    ok scalar(@glim_exons), $num_exons[$i];
    if($i == 5) {
	ok $glim_exons[1]->start, 30152;
	ok $glim_exons[1]->end, 30235;
	ok $glim_exons[1]->strand, -1;
    }
}

# Glimmer testing (GlimmerM)
my $ghmm = Bio::Tools::Glimmer->new('-file' => Bio::Root::IO->catfile(qw(t data GlimmerHMM.out)));
my $ghmmgene = $ghmm->next_prediction;

ok($ghmmgene);
ok($ghmmgene->seq_id, 'gi|23613028|ref|NC_004326.1|');
ok($ghmmgene->source_tag, 'GlimmerHMM');
ok($ghmmgene->primary_tag, 'transcript');
ok($ghmmgene->exons == 1);

@num_exons = qw(0 1 2 4 2 2 1 1 1 2 2 2 10 4 1 1); # only first few tested
$i = 1;
while ($ghmmgene = $ghmm->next_prediction) {
  $i++;
  my @ghmm_exons = $ghmmgene->exons;    
  ok(scalar(@ghmm_exons), $num_exons[$i]) if $i <= $#num_exons;
  if ($i == 9) {
    ok( $ghmm_exons[1]->start, 5538 );
    ok( $ghmm_exons[1]->end,   5647 );
    ok( $ghmm_exons[1]->strand > 0  );
  }
}
ok($i, 44);


