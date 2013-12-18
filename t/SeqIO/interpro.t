# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  
  test_begin(-tests => 20,
             -requires_module => 'XML::DOM::XPath');
  
  use_ok('Bio::SeqIO::interpro');
}

my $verbose = test_debug();

my $t_file = test_input_file('test.interpro');
my $a_in = Bio::SeqIO->new( -file    => $t_file,
                            -verbose => $verbose,
                            -format  => 'interpro');
isa_ok($a_in, 'Bio::SeqIO');

my $seq = $a_in->next_seq();
ok($seq, 'seq obj is defined');
isa_ok($seq, 'Bio::Seq::RichSeq');
is(scalar( $seq->get_SeqFeatures() ), 6, 'right number of SeqFeatures');

my($feat) = $seq->get_SeqFeatures();
isa_ok($feat,'Bio::SeqFeature::Generic');

is($feat->display_name,'Retinoblastoma-associated protein, B-box', 'display_name()');

ok($seq = $a_in->next_seq(), 'seq object is defined');
is(scalar( $seq->get_SeqFeatures() ),40, 'right number of SeqFeatures');

ok(!($seq = $a_in->next_seq()), 'there is no next_seq (correctly)');

# Bug 1908 (enhancement)
$t_file = test_input_file('interpro_ebi.xml');
my $b_in = Bio::SeqIO->new( -file    => $t_file,
                            -verbose => $verbose,
                            -format  => 'interpro');
$seq = $b_in->next_seq();
ok($seq, 'bug 1908');

my @features = $seq->get_SeqFeatures;
is(scalar @features, 2, 'right number of SeqFeatures');
is($features[0]->primary_tag, 'region', 'primary_tag()');
is($features[0]->display_name, 'Protein of unknown function DUF1021', 'display_name()');
is($features[0]->location->end, 78, 'location->end()');

my @dblinks = $features[0]->annotation->get_Annotations('dblink');
is(scalar @dblinks, 3, 'right number of dblinks');
is($dblinks[1]->primary_id, 'IPR009366', 'first primary_id');
is($dblinks[2]->primary_id, 'PF06257.1', 'second primary_id');

my $other_t_file = test_input_file('test.interpro-go.xml');
my $ipr_in = Bio::SeqIO->new( -file    => $other_t_file,
                              -verbose => $verbose,
                              -format  => 'interpro');

$seq = $ipr_in->next_seq();
@features = $seq->get_SeqFeatures;
@dblinks = $features[0]->annotation->get_Annotations('dblink');
is(scalar @dblinks, 4, 'right number of dblinks');
is($dblinks[3]->primary_id, 'GO:0003677', 'primary_id via dblinks');
