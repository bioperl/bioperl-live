# -*-Perl-*- Test Harness script for Bioperl
# $Id: arp.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests           => 48,
               -requires_module => 'Data::Stag');
	
	use_ok('Bio::AlignIO::arp');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# ARP format
$str  = Bio::AlignIO ->new(
    '-file'	=> test_input_file("testaln.arp"),
    -verbose => 1,
    '-format'	=> 'arp');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '01/1-399','ARP get_nse()');
is($aln->get_seq_by_pos(1)->length, '407');
is($aln->num_sequences, 60,'ARP num_sequences()');
is($aln->id, 'Mandenka', 'ARP id()');
is($aln->description, 'mtDNA sequences in the Senegalese Mandenka (hypervariable region 1)', 'ARP description()');
my $coll = $aln->annotation;
isa_ok($coll, 'Bio::AnnotationCollectionI');
my ($ann) = $coll->get_Annotations('Samples');
isa_ok($ann, 'Bio::AnnotationI');
my %nodes = $ann->pairs;
is(keys %nodes, 60);
is($nodes{'03'}, 10);
is(($coll->get_Annotations('DataType'))[0]->value,'DNA');
is(($coll->get_Annotations('MissingData'))[0]->value,'?');

$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("testaln2.arp"),
    '-format'	=> 'arp');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '000/1-29','ARP get_nse()');
is($aln->num_sequences, 3,'ARP num_sequences()');
is($aln->id, 'Population 1', 'ARP id()');
is($aln->description, 'An example of DNA sequence data', 'ARP description()');
$coll = $aln->annotation;
isa_ok($coll, 'Bio::AnnotationCollectionI');
($ann) = $coll->get_Annotations('Samples');
isa_ok($ann, 'Bio::AnnotationI');
%nodes = $ann->pairs;
is(keys %nodes, 3);
is($nodes{'001'}, 1);
is(($coll->get_Annotations('DataType'))[0]->value, 'DNA');
is(($coll->get_Annotations('SampleSize'))[0]->value, 6);

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->get_nse, '001/1-29','ARP get_nse()');
is($aln->num_sequences, 8,'ARP num_sequences()');
is($aln->id, 'Population 2', 'ARP id()');
is($aln->description, 'An example of DNA sequence data', 'ARP description()');
$coll = $aln->annotation;
isa_ok($coll, 'Bio::AnnotationCollectionI');
($ann) = $coll->get_Annotations('Samples');
isa_ok($ann, 'Bio::AnnotationI');
%nodes = $ann->pairs;
is(keys %nodes, 8);
is($nodes{'001'}, 1);
is(($coll->get_Annotations('DataType'))[0]->value, 'DNA');
is(($coll->get_Annotations('SampleSize'))[0]->value, 8);

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->get_nse, '024/1-29','ARP get_nse()');
is($aln->num_sequences, 6,'ARP num_sequences()');
is($aln->id, 'Population 3', 'ARP id()');
is($aln->description, 'An example of DNA sequence data', 'ARP description()');
$coll = $aln->annotation;
isa_ok($coll, 'Bio::AnnotationCollectionI');
($ann) = $coll->get_Annotations('Samples');
isa_ok($ann, 'Bio::AnnotationI');
%nodes = $ann->pairs;
is(keys %nodes, 6);
is($nodes{'024'}, 1);
is(($coll->get_Annotations('DataType'))[0]->value, 'DNA');
is(($coll->get_Annotations('SampleSize'))[0]->value, 6);
