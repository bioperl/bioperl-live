# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 28);

    use_ok('Bio::SeqIO');
    use_ok('Bio::SeqFeature::Gene::Transcript');
    use_ok('Bio::SeqFeature::Gene::UTR');
    use_ok('Bio::SeqFeature::Gene::Exon');
    use_ok('Bio::SeqFeature::Gene::Poly_A_site');
    use_ok('Bio::SeqFeature::Gene::GeneStructure');
    use_ok('Bio::Location::Fuzzy');
}


my ( $seqio, $geneseq, $gene, $transcript, $poly_A_site1, $poly_A_site2,
     $fiveprimeUTR, $exon);

# tests for Bio::SeqFeature::Gene::* objects
# using information from acc: AB077698 as a guide

ok $seqio = Bio::SeqIO->new(
    -format => 'genbank',
    -file   => test_input_file('AB077698.gb'),
);

ok $geneseq = $seqio->next_seq();

ok $gene = Bio::SeqFeature::Gene::GeneStructure->new(
    -primary => 'gene',
    -start   => 1,
    -end     => 2701,
    -strand  => 1,
);

ok $transcript = Bio::SeqFeature::Gene::Transcript->new(
    -primary => 'CDS',
    -start   => 80,
    -end     => 1144,
    -tag     => { 
        'gene' => "CHCR",
        'note' => "Cys3His CCG1-Required Encoded on BAC clone RP5-842K24 (AL050310) The human CHCR (Cys3His CCG1-Required) protein is highly related to EXP/MBNL (Y13829, NM_021038, AF401998) and MBLL (NM_005757,AF061261), which together comprise the human Muscleblind family",
       'codon_start' => 1,
       'protein_id'  => 'BAB85648.1',
    }
);

ok $poly_A_site1 = Bio::SeqFeature::Gene::Poly_A_site->new(
    -primary => 'polyA_site',
    -start => 2660,
    -end   => 2660,
    -tag   => { 
        'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#2 used by CHCR EST clone DKFZp434G2222 (AL133625)"
    }
);

ok $poly_A_site2 = Bio::SeqFeature::Gene::Poly_A_site->new(
    -primary => 'polyA_site',
    -start => 1606,
    -end   => 1606,
    -tag   => { 
        'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#1 used by CHCR EST clone PLACE1010202 (AK002178)",
    }
);

ok $fiveprimeUTR = Bio::SeqFeature::Gene::UTR->new(-primary => "utr5prime");
ok $fiveprimeUTR->location(
    Bio::Location::Fuzzy->new(
        -start => "<1",
        -end   => 79,
    )
);
ok my $threeprimeUTR = Bio::SeqFeature::Gene::UTR->new(
    -primary => "utr3prime",
    -start   => 1145,
    -end     => 2659,
);

# Did a quick est2genome against genomic DNA (this is on Chr X) to
# get the gene structure by hand since it is not in the file
# --Jason

ok $exon = Bio::SeqFeature::Gene::Exon->new(
    -primary => 'exon',
    -start => 80,
    -end   => 177,
);
ok $geneseq->add_SeqFeature($exon);

ok $geneseq->add_SeqFeature($fiveprimeUTR);
ok $geneseq->add_SeqFeature($threeprimeUTR);
ok $geneseq->add_SeqFeature($poly_A_site1);
ok $geneseq->add_SeqFeature($poly_A_site2);

ok $transcript->add_utr($fiveprimeUTR, 'utr5prime');
ok $transcript->add_utr($threeprimeUTR, 'utr3prime');

ok $transcript->add_exon($exon);

# API only supports a single poly-A site per transcript at this point 
$transcript->poly_A_site($poly_A_site2);
$geneseq->add_SeqFeature($transcript);
$gene->add_transcript($transcript);
$geneseq->add_SeqFeature($gene);

my ($t) = $gene->transcripts(); # get 1st transcript
ok(defined $t); 
is($t->mrna->length, 1693, 'mRNA spliced length');
is($gene->utrs, 2, 'has 2 UTRs');


