# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 10);
	
	use_ok('Bio::AlignIO::fasta');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# FASTA
$str = Bio::AlignIO->new(
		 -file => test_input_file("testaln.fasta"), 
		 -format => 'fasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431/1-318', 
  "fasta input test ";
is ($aln->get_seq_by_pos(1)->description, 'DESCRIPTION HERE', 
    "fasta input test for description");
is ($aln->get_seq_by_pos(11)->display_id, 'AK_YEAST',
    "fasta input test for id");

is ($aln->get_seq_by_pos(2)->end, 318,
    "fasta input test for end");

is ($aln->get_seq_by_pos(11)->description, 'A COMMENT FOR YEAST', 
    "fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'fasta');
$status = $strout->write_aln($aln);
is $status, 1,"fasta output test";

my $in = Bio::AlignIO->newFh(
   '-file'  => test_input_file("testaln.fasta"), 
			       '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(
   '-file' => ">".test_output_file(), 
				'-format' => 'pfam');
while ( $aln = <$in>) {
    is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431/1-318',
     "filehandle input test  ";
    $status = print $out $aln;
    last;
}
is $status, 1, "filehandle output test";
