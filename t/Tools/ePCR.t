# -*-Perl-*- Test Harness script for Bioperl
# $Id$


BEGIN {     
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 27);
	
    use_ok('Bio::Tools::EPCR');
    use_ok('Bio::SeqIO');
}

my $DEBUG = test_debug();

my $seqio = Bio::SeqIO->new('-format' => 'fasta', '-file' => test_input_file('genomic-seq.fasta'));

my $seq = $seqio->next_seq;
ok($seq);
my $epcr = Bio::Tools::EPCR->new( '-file' => test_input_file('genomic-seq.epcr'));
ok ($epcr);
my %strand;
while( defined(my $feature = $epcr->next_feature) ) {
    ok($feature);
    ok($feature->start);
    ok($feature->end);
    $seq->add_SeqFeature($feature);
    $strand{$feature->strand} ++;
}
is ($strand{1},  3, 'got 3 forward strand ePCR hits');
is ($strand{-1}, 3, 'got 3 reverse strand ePCR hits');

if( $DEBUG ) {
    $seqio = Bio::SeqIO->new('-format' => 'genbank' );
    $seqio->write_seq($seq);
}
