# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
  use lib 't/lib';
  use BioperlTest;
  
  test_begin(-tests => 17,
			 -requires_module => 'XML::DOM::XPath');
  
  use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

my $t_file = test_input_file('test.interpro');
my $a_in = Bio::SeqIO->new( -file => $t_file,
									 -verbose => $verbose,
									 -format => 'interpro');

my $seq = $a_in->next_seq();
ok($seq);
isa_ok($seq, 'Bio::Seq::RichSeq');
is(scalar( $seq->get_SeqFeatures() ),6);

my($feat) = $seq->get_SeqFeatures();
isa_ok($feat,'Bio::SeqFeature::Generic');

is($feat->display_name,'Retinoblastoma-associated protein, B-box');

ok($seq = $a_in->next_seq());
is(scalar( $seq->get_SeqFeatures() ),40);

ok(!($seq = $a_in->next_seq()));

# Bug 1908 (enhancement)
$t_file = test_input_file('interpro_ebi.xml');
my $b_in = Bio::SeqIO->new( -file => $t_file,
									 -verbose => $verbose,
									 -format => 'interpro');
$seq = $b_in->next_seq();
ok($seq);

my @features = $seq->get_SeqFeatures;
is scalar @features,2;
is $features[0]->primary_tag, 'region';
is $features[0]->display_name,'Protein of unknown function DUF1021';
is $features[0]->location->end,78;

my @dblinks = $features[0]->annotation->get_Annotations('dblink');
is (scalar @dblinks,3);
is $dblinks[1]->primary_id,'IPR009366';
is $dblinks[2]->primary_id,'PF06257.1';
