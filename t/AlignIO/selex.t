# -*-Perl-*- Test Harness script for Bioperl
# $Id: selex.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 4);
	
	use_ok('Bio::AlignIO::selex');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# SELEX
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.selex"),
			   '-format' => 'selex');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'HSFAU/1-518', "selex format test ";

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'selex');
$status = $strout->write_aln($aln);
is $status, 1, "selex output test";
