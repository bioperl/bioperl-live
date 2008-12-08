# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 185);

    use_ok('Bio::Tools::Fgenesh');
    use_ok('Bio::Tools::Genscan');
    use_ok('Bio::Tools::Genemark');
    use_ok('Bio::Tools::Glimmer');
    use_ok('Bio::Tools::MZEF');  
    use_ok('Bio::SeqIO');
}

# Genscan report
my $genscan = Bio::Tools::Genscan->new('-file' => test_input_file('genomic-seq.genscan'));
ok $genscan;

# original sequence
my $seqin = Bio::SeqIO->new('-file' => test_input_file('genomic-seq.fasta'),
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
	is $fea->strand(), -1, 
	     "strand match (".$fea->strand()."  and -1)";
	$fea = ($gene->poly_A_site());
	is $fea->score(), 1.05, 
             "score match (".$fea->score()." and 1.05)";
    }
    if($pred_num == 2) {
	$fea = ($gene->exons("Initial"))[0];
	is $fea->strand(), 1, 
	"strand match (".$fea->strand()." and 1)";
	is $fea->score(), 4.46, 
             "score match (".$fea->score()." and 4.46)";
    }
    if($pred_num == 3) {
	my @exons = $gene->exons("Initial");
	is scalar(@exons), 0, 
	     "initial exons ".scalar(@exons);
	$fea = ($gene->exons())[0];
	is $fea->score(),  1.74, 
             "score match ".$fea->score();
    }
    if($seq) {
	$prtseq = $gene->predicted_protein()->seq();
        $cds = $gene->cds();
	ok($cds);
	$tr_cds = $cds->translate()->seq();
	$tr_cds =~ s/\*$//;
	is( lc($prtseq), lc($tr_cds),
	    "predicted and extracted protein seqs match");
    }
}

# Genscan report with no genes predicted
my $null_genscan = Bio::Tools::Genscan->new('-file' => test_input_file('no-genes.genscan'));
ok $null_genscan;
my $no_gene = $null_genscan->next_prediction;
my @exons = $no_gene->exons;
is($#exons,-1);

# MZEF report
my $mzef = Bio::Tools::MZEF->new('-file' => test_input_file('genomic-seq.mzef'));
ok $mzef;

my $exon_num = 0;
my $gene = $mzef->next_prediction();

is($gene->exons, 23);

# Genemark testing
my $genemark = Bio::Tools::Genemark->new('-file' => test_input_file('genemark.out'));

my $gmgene = $genemark->next_prediction();
is $gmgene->seq_id(), "Hvrn.contig8";
is $genemark->analysis_date(), "Thu Mar 22 10:25:00 2001";

my $i = 0;
my @num_exons = (1,5,2,1,9,5,3,2,3,2,1,2,7);
while($gmgene = $genemark->next_prediction()) {
    $i++;
    my @gmexons = $gmgene->exons();
    is scalar(@gmexons), $num_exons[$i];

    if($i == 5) {
	my $gmstart = $gmexons[0]->start();
	is $gmstart, 23000;

	my $gmend = $gmexons[0]->end();
	is $gmend, 23061;
    }
}

# Genemark testing (prokaryotic gene fragment)
$genemark = Bio::Tools::Genemark->new('-file'    => test_input_file('genemark-fragment.out'),
                                         '-seqname' => 'AAVN02000021.1');

$gmgene = $genemark->next_prediction();
is $gmgene->seq_id(), 'AAVN02000021.1','Genemark tests';
is $gmgene->start(), 2;
is $gmgene->end(), 214;
is $gmgene->strand(), '1';
my ($gmexon) = $gmgene->exons();
isa_ok $gmexon->location(), 'Bio::Location::Fuzzy';
is $gmexon->location->start_pos_type(), 'BEFORE';
is $gmexon->location->max_start(), 2;
is $gmexon->location->end_pos_type(), 'EXACT';
is $gmexon->location->end(), 214;

$gmgene = $genemark->next_prediction();
is $gmgene->seq_id(), 'AAVN02000021.1';
is $gmgene->start, 459;
is $gmgene->end, 596;
is $gmgene->strand(), '1';
($gmexon) = $gmgene->exons();
isa_ok $gmexon->location, 'Bio::Location::Fuzzy';
is $gmexon->location->start_pos_type(), 'EXACT';
is $gmexon->location->start(), 459;
is $gmexon->location->end_pos_type(), 'AFTER';
is $gmexon->location->min_end(), 596;

# Glimmer testing (GlimmerM)
my $glimmer_m = Bio::Tools::Glimmer->new('-file' => test_input_file('GlimmerM.out'));
$gmgene = $glimmer_m->next_prediction;

ok($gmgene);
is($gmgene->seq_id, 'gi|23613028|ref|NC_004326.1|');
is($gmgene->source_tag, 'GlimmerM_3.0');
is($gmgene->primary_tag, 'transcript');
is(($gmgene->get_tag_values('Group'))[0], 'GenePrediction1');
my @glim_exons = $gmgene->exons;
is(scalar (@glim_exons), 1);
is($glim_exons[0]->start, 461);
is($glim_exons[0]->end, 523);
is($glim_exons[0]->strand, -1);
is(($glim_exons[0]->get_tag_values('Group'))[0], 'GenePrediction1');

@num_exons = (0,1,3,1,4,2,5,2,8,3,5);
$i = 1;
while($gmgene = $glimmer_m->next_prediction()) {
    $i++;
    is(($gmgene->get_tag_values('Group'))[0],"GenePrediction$i");
    @glim_exons = $gmgene->exons();    
    is scalar(@glim_exons), $num_exons[$i];
    if($i == 5) {
	is $glim_exons[1]->start, 23910;
	is $glim_exons[1]->end, 23956;
	is $glim_exons[1]->strand, 1;
    }
}

# Glimmer testing (GlimmerHMM)
my $glimmer_hmm = Bio::Tools::Glimmer->new('-file' => test_input_file('GlimmerHMM.out'));
my $ghmmgene = $glimmer_hmm->next_prediction;

ok($ghmmgene);
is($ghmmgene->seq_id, 'gi|23613028|ref|NC_004326.1|');
is($ghmmgene->source_tag, 'GlimmerHMM');
is($ghmmgene->primary_tag, 'transcript');
is($ghmmgene->exons, 1);

@num_exons = qw(0 1 2 4 2 2 1 1 1 2 2 2 10 4 1 1); # only first few tested
$i = 1;
while ($ghmmgene = $glimmer_hmm->next_prediction) {
  $i++;
  my @ghmm_exons = $ghmmgene->exons;    
  is(scalar(@ghmm_exons), $num_exons[$i]) if $i <= $#num_exons;
  if ($i == 9) {
    is( $ghmm_exons[1]->start, 5538 );
    is( $ghmm_exons[1]->end,   5647 );
    cmp_ok( $ghmm_exons[1]->strand, '>', 0  );
  }
}
is($i, 44);

# Glimmer testing (Glimmer 2.X)
my $glimmer_2 = Bio::Tools::Glimmer->new('-file' => test_input_file('Glimmer2.out'),
										 '-seqname' => 'BCTDNA',
										 '-seqlength' => 29940,);
my $g2gene = $glimmer_2->next_prediction;

ok($g2gene);
is($g2gene->seq_id, 'BCTDNA');
is($g2gene->source_tag, 'Glimmer_2.X');
is($g2gene->primary_tag, 'gene');
is($g2gene->start, 292);
is($g2gene->end, 1623);
is($g2gene->frame, 0);
is($g2gene->strand, 1);

$i = 1;
while ($g2gene = $glimmer_2->next_prediction) {
    $i++;
    if ($i == 2) {
        is($g2gene->start, 2230);
        is($g2gene->end, 2349);
        is($g2gene->strand, -1);
        is($g2gene->frame, 0);		
	} elsif ($i == 25) {
        isa_ok($g2gene->location, 'Bio::Location::SplitLocationI');
        my @sublocations = $g2gene->location->sub_Location();
        is(scalar (@sublocations), 2);
        is($sublocations[0]->start, 29263);
        is($sublocations[0]->end, 29940);
        is($sublocations[1]->start, 1);
        is($sublocations[1]->end, 9);
        is($g2gene->strand, 1);
        is($g2gene->frame, 0);
    }
}
is($i, 25);


# Glimmer testing (Glimmer 3.X)
my $glimmer_3 = Bio::Tools::Glimmer->new('-file' => test_input_file('Glimmer3.predict'),
										 '-detail' => test_input_file('Glimmer3.detail'));
my $g3gene = $glimmer_3->next_prediction;

ok($g3gene);
is($g3gene->seq_id, 'BCTDNA');
is($g3gene->source_tag, 'Glimmer_3.X');
is($g3gene->primary_tag, 'gene');
is($g3gene->score, '9.60');
isa_ok($g3gene->location, 'Bio::Location::SplitLocationI');
my @sublocations = $g3gene->location->sub_Location();
is(scalar (@sublocations), 2);
is($sublocations[0]->start, 29263);
is($sublocations[0]->end, 29940);
is($sublocations[1]->start, 1);
is($sublocations[1]->end, 9);
is($g3gene->frame, 0);

$i = 1;
while ($g3gene = $glimmer_3->next_prediction) {
    $i++;
    if ($i == 13) {
        is($g3gene->start, 13804);
        is($g3gene->end, 14781);
        is($g3gene->strand, -1);
        is($g3gene->frame, 0);
        is($g3gene->score, '5.51');
		
        my ($orfid) = $g3gene->has_tag('Group') ? $g3gene->get_tag_values('Group') : undef;
        is($orfid, 'GenePrediction_00015');
    }
}
is($i, 27);

# Glimmer 3.X (prokaryotic gene fragment)
my $glimmer_3a = Bio::Tools::Glimmer->new(
                                         '-file'   => test_input_file('glimmer3-fragment.predict'),
                                         '-detail' => test_input_file('glimmer3-fragment.detail'), 
                                        );
my $g3gene_a = $glimmer_3a->next_prediction;

ok($g3gene_a);

isa_ok $g3gene_a->location(), 'Bio::Location::Fuzzy';
is $g3gene_a->location->start_pos_type(), 'BEFORE';
is $g3gene_a->location->max_start(), 1;
is $g3gene_a->location->end_pos_type(), 'EXACT';
is $g3gene_a->location->end(), 674;
is $g3gene_a->frame(), 2;

for (1..3) { $g3gene_a = $glimmer_3a->next_prediction; }

isa_ok $g3gene_a->location(), 'Bio::Location::Fuzzy';
is $g3gene_a->location->start_pos_type(), 'EXACT';
is $g3gene_a->location->start(), 2677;
is $g3gene_a->frame(), 0;
is $g3gene_a->location->end_pos_type(), 'AFTER';
is $g3gene_a->location->min_end(), 2932;
is $g3gene_a->score, '5.63';

# Fgenesh
my $fgh = Bio::Tools::Fgenesh->new(
                                   '-file' => test_input_file('fgenesh.out'),
                                  );

my $fghgene = $fgh->next_prediction();

ok($fghgene);
is($fghgene->seq_id, 'gi|1914348|emb|Z81551.1|');
is($fghgene->source_tag, 'Fgenesh');
is($fghgene->start(), 29);
is($fghgene->end(), 1869);
cmp_ok($fghgene->strand(), '<', 0);

$i = 0;
@num_exons = (2,5,4,8);

while ($fghgene = $fgh->next_prediction()) {

    $i++;
    my @fghexons = $fghgene->exons();
    is(scalar(@fghexons), $num_exons[$i]);

    if ($i == 2) {
        cmp_ok($fghexons[0]->strand(), '>', 0);
        is($fghexons[0]->primary_tag(), 'InitialExon');
        is($fghexons[0]->start(), 14778);
        is($fghexons[0]->end(), 15104);
        cmp_ok($fghexons[3]->strand(), '>', 0);
        is($fghexons[3]->primary_tag(), 'TerminalExon');        
        is($fghexons[3]->start(), 16988);
        is($fghexons[3]->end(), 17212);
    }

}

is($i, 3);
