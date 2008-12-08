# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
	use_ok('Bio::Factory::SeqAnalysisParserFactory');
	use_ok('Bio::SeqIO');
}

my ($seqio,$seq,$factory,$parser, $gene_seen, $exon_seen);

$seqio = Bio::SeqIO->new('-format'=>'fasta', '-file' => test_input_file('genomic-seq.fasta'));
isa_ok $seqio, 'Bio::SeqIO';
$seq = $seqio->next_seq;
isa_ok $seq, 'Bio::PrimarySeqI';

$factory = Bio::Factory::SeqAnalysisParserFactory->new();

# let's test the genscan factory
$parser = $factory->get_parser(-input => test_input_file('genomic-seq.genscan'),
				  -method => 'genscan');
isa_ok $parser, 'Bio::SeqAnalysisParserI';
while( my $feat = $parser->next_feature() ){
    $seq->add_SeqFeature($feat);
}
($gene_seen, $exon_seen)  = (0,0);
foreach my $feat (  $seq->top_SeqFeatures() ) {
    if( $feat->isa("Bio::Tools::Prediction::Gene") ) {
	foreach my $exon ( $feat->exons ) {
	    $exon_seen++;
	}
	$gene_seen++;
    } 
}
is $exon_seen, 37;
is $gene_seen, 3;

# let's test the mzef factory
$parser = $factory->get_parser(-input => test_input_file('genomic-seq.mzef'),
			       -method=> 'mzef');
$seqio = Bio::SeqIO->new('-format'=>'fasta', '-file' => test_input_file('genomic-seq.fasta'));
ok $seq = $seqio->next_seq();
isa_ok $seq, 'Bio::PrimarySeqI';

isa_ok $parser, 'Bio::SeqAnalysisParserI';
while( my $feat = $parser->next_feature() ){
    $seq->add_SeqFeature($feat);
}
($gene_seen, $exon_seen)  = (0,0);
foreach my $feat (  $seq->top_SeqFeatures() ) {
    if( $feat->isa("Bio::Tools::Prediction::Gene") ) {
	foreach my $exon ( $feat->exons ) { 
	    $exon_seen++;
	}
	$gene_seen++;
    } 
}
is $exon_seen, 23;
is $gene_seen, 1;

# let's test the ePCR factory

$parser = $factory->get_parser(-input => test_input_file('genomic-seq.epcr'),
			       -method => 'epcr');

$seq->flush_SeqFeatures;

isa_ok $parser, 'Bio::SeqAnalysisParserI';
while( my $feat = $parser->next_feature() ){
    $seq->add_SeqFeature($feat);
}

is $seq->top_SeqFeatures(), 7;
