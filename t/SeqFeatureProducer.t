# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$

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
    plan tests => 9;
}

use Bio::SeqFeatureProducer;
use Bio::Tools::MZEF;
use Bio::SeqIO;

my ($seqio,$seq,$sfp, $gene_seen, $exon_seen);

$seqio = new Bio::SeqIO('-format'=>'fasta', '-file' => 't/genomic-seq.fasta');
ok $seqio->isa('Bio::SeqIO');# 'seqio was not created';
$seq = $seqio->next_seq;
ok $seq->isa('Bio::PrimarySeqI');#'could not read sequence';

$sfp = new Bio::SeqFeatureProducer(-method => 'genscan',
					     -input  => 't/genomic-seq.genscan');
ok $sfp->isa('Bio::SeqFeatureProducer');#'no SeqFeatureProducer created';

$sfp->add_features($seq);
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
 
$sfp = new Bio::SeqFeatureProducer();
ok($sfp->isa('Bio::SeqFeatureProducer'));

my $parser = new Bio::Tools::MZEF(-file => 't/genomic-seq.mzef');
$seqio = new Bio::SeqIO('-format'=>'fasta', '-file' => 't/genomic-seq.fasta');

$seq = $seqio->next_seq();
ok(defined $seq && $seq->isa('Bio::PrimarySeqI'));

$sfp->add_features($seq,$parser);

($gene_seen, $exon_seen)  = (0,0);
foreach my $feat ( $seq->top_SeqFeatures() ) {
    if( $feat->isa("Bio::Tools::Prediction::Gene") ) {
	foreach my $exon ( $feat->exons ) {
	    $exon_seen++;
	}
	$gene_seen++;
    } 
}
ok $exon_seen, 23;
ok $gene_seen, 1;
