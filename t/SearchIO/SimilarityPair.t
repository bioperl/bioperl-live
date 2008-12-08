# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 12);
	
	use_ok('Bio::SearchIO');
	use_ok('Bio::SeqIO');
}

# test SimilarityPair

ok my $seq = (Bio::SeqIO->new('-format' => 'fasta',
			  '-file' => test_input_file('AAC12660.fa')))->next_seq();
isa_ok $seq, 'Bio::SeqI';
ok my $blast = Bio::SearchIO->new('-file'=>test_input_file('blast.report'), '-format' => 'blast');
isa_ok $blast, 'Bio::SearchIO';
my $r = $blast->next_result;
ok my $hit = $r->next_hit;
isa_ok $hit, 'Bio::Search::Hit::HitI';
ok my $sim_pair = $hit->next_hsp;
isa_ok $sim_pair, 'Bio::SeqFeatureI';
ok $seq->add_SeqFeature($sim_pair);
is $seq->all_SeqFeatures(), 1;
