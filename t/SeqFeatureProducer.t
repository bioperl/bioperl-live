# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


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
BEGIN { $| = 1; print "1..10\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::SeqFeatureProducer;
use Bio::Tools::MZEF;
use Bio::SeqIO;
$loaded = 1;
print "ok 1\n";    # 1st test passes.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

my ($seqio,$seq,$sfp, $gene_seen, $exon_seen);

test 2, defined ( $seqio = new Bio::SeqIO(-format=>'fasta', -file => 't/genomic-seq.fasta')), 'seqio was not created';
test 3, defined ( $seq = $seqio->next_seq ), 'could not read sequence';

test 4, ( $sfp = new Bio::SeqFeatureProducer(-method => 'genscan',
					     -input  => 't/genomic-seq.genscan')), 'no SeqFeatureProducer created';

$sfp->add_features($seq);
($gene_seen, $exon_seen)  = (0,0);
foreach my $feat ( @features = $seq->top_SeqFeatures() ) {
    if( $feat->isa("Bio::Tools::Prediction::Gene") ) {
	foreach my $exon ( $feat->exons ) {
	    $exon_seen++;
	}
	$gene_seen++;
    } 
}
test 5, $exon_seen == 37;
test 6, $gene_seen == 3;

;
test 7, defined ( $sfp = new Bio::SeqFeatureProducer());

my $parser = new Bio::Tools::MZEF(-file => 't/genomic-seq.mzef');
$seqio = new Bio::SeqIO(-format=>'fasta', -file => 't/genomic-seq.fasta');

test 8, defined ($seq = $seqio->next_seq());

$sfp->add_features($seq,$parser);

($gene_seen, $exon_seen)  = (0,0);
foreach my $feat ( @features = $seq->top_SeqFeatures() ) {
    if( $feat->isa("Bio::Tools::Prediction::Gene") ) {
	foreach my $exon ( $feat->exons ) {
	    $exon_seen++;
	}
	$gene_seen++;
    } 
}
test 9, $exon_seen == 23;
test 10, $gene_seen == 1;
