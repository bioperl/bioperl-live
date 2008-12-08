# -*-Perl-*- Test Harness script for Bioperl
# $Id: pfam.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::AlignIO::pfam');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PFAM format (no annotation)
$str = Bio::AlignIO->new(
	  '-file' => test_input_file("testaln.pfam"));
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			    '-format' => 'pfam');
$status = $strout->write_aln($aln);
is($status, 1, " pfam output test");

