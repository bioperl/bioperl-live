# -*-Perl-*-

use strict;
BEGIN { 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 11;
}

use Bio::Factory::SeqAnalysisParserFactory;
use Bio::SeqIO;
use Bio::Root::IO;

my ($seqio,$seq,$factory,$parser, $gene_seen, $exon_seen);

$seqio = new Bio::SeqIO('-format'=>'fasta', '-file' => Bio::Root::IO->catfile("t","data","genomic-seq.fasta"));
ok $seqio->isa('Bio::SeqIO');# 'seqio was not created';
$seq = $seqio->next_seq;
ok $seq->isa('Bio::PrimarySeqI');#'could not read sequence';

$factory = new Bio::Factory::SeqAnalysisParserFactory();

# let's test the genscan factory
$parser = $factory->get_parser(-input => Bio::Root::IO->catfile("t","data","genomic-seq.genscan"),
				  -method => 'genscan');
ok $parser->isa('Bio::SeqAnalysisParserI');#'noSeqAnalysisParserI created';
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
ok $exon_seen, 37;
ok $gene_seen, 3;

# let's test the mzef factory
$parser = $factory->get_parser(-input => Bio::Root::IO->catfile("t","data","genomic-seq.mzef"),
			       -method=> 'mzef');
$seqio = new Bio::SeqIO('-format'=>'fasta', '-file' => Bio::Root::IO->catfile("t","data","genomic-seq.fasta"));
$seq = $seqio->next_seq();
ok(defined $seq && $seq->isa('Bio::PrimarySeqI'));

ok $parser->isa('Bio::SeqAnalysisParserI');#'noSeqAnalysisParserI created';
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
ok $exon_seen, 23;
ok $gene_seen, 1;

# let's test the ePCR factory

$parser = $factory->get_parser(-input => Bio::Root::IO->catfile("t","data", "genomic-seq.epcr"),
			       -method => 'epcr');

$seq->flush_SeqFeatures;

ok $parser->isa('Bio::SeqAnalysisParserI');#'noSeqAnalysisParserI created';
while( my $feat = $parser->next_feature() ){
    $seq->add_SeqFeature($feat);
}

ok $seq->top_SeqFeatures(), 7;
