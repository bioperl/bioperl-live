# -*-Perl-*- Test Harness script for Bioperl
# $Id: nexus.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 43);
	
	use_ok('Bio::AlignIO::nexus');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# NEXUS
$str = Bio::AlignIO->new(
   '-file' => test_input_file('testaln.nexus'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapiens/1-45';
$strout = Bio::AlignIO->new('-file'  => ">".test_output_file(),
			    '-format' => 'nexus', );
$status = $strout->write_aln($aln);
is $status, 1, "nexus output test";

$str = Bio::AlignIO->new(
   '-file' => test_input_file('Bird_Ovomucoids.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('basic-ladder.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('Kingdoms_DNA.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file('char-interleave.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('Primate_mtDNA.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("char-matrix-spaces.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family4nl.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family7n.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("long-names.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family8a.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("multiline-intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("Treebase-chlamy-dna.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("quoted-strings1.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("UnaSmithHIV-both.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("quoted-strings2.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("barns-combined.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("radical-whitespace.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("basic-bush.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("radical-whitespace_02.nex"),
			  '-format' => 'nexus');
