# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 13);
	
	use_ok('Bio::Tools::EMBOSS::Palindrome');
	use_ok('Bio::Tools::GFF');
}

my $DEBUG = test_debug();

my $parser = Bio::Tools::EMBOSS::Palindrome->new(-file => test_input_file('humts1.pal'));

my $seq = $parser->next_seq;
ok($seq);
is($seq->display_id, 'HUMTS1');
is($seq->length, 18596);
my @features = $seq->get_SeqFeatures();
is(scalar @features, 23);

is($features[0]->feature1->start, 126);
is($features[0]->feature1->end, 142);
is($features[0]->feature1->strand, 1);
is($features[0]->feature1->seq_id, 'HUMTS1');


is($features[0]->feature2->start, 201);
is($features[0]->feature2->end, 217);
is($features[0]->feature2->strand, -1);

if( $DEBUG ) {
    my $out = Bio::Tools::GFF->new(-gff_version => 2);
    $out->write_feature($features[0]);
}
