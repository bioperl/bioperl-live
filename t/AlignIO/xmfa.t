# -*-Perl-*- Test Harness script for Bioperl
# $Id: xmfa.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
	use_ok('Bio::AlignIO::xmfa');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# XMFA
$str = Bio::AlignIO->new(
		 -file => test_input_file("testaln.xmfa"), 
		 -format => 'xmfa');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'chrY/1-598', 
  "xmfa input test ";
is $aln->get_seq_by_pos(1)->strand, 1, 
  "xmfa strand test";
is ($aln->get_seq_by_pos(2)->description, undef, 
    "xmfa input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'chr7',
    "xmfa input test for id");
is ($aln->get_seq_by_pos(2)->start, 5000,
    "xmfa input test for end");
is ($aln->get_seq_by_pos(2)->end, 5534,
    "xmfa input test for end");
is ($aln->score, 111, 'xmfa alignment score');

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'chrY/1000-1059', 
  "xmfa input test ";
is $aln->get_seq_by_pos(1)->strand, 1, 
  "xmfa strand";
is ($aln->get_seq_by_pos(2)->description, undef, 
    "xmfa input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'chr12',
    "xmfa input test for id");
is ($aln->get_seq_by_pos(2)->start, 6000,
    "xmfa input test for end");
is ($aln->get_seq_by_pos(1)->end, 1059,
    "xmfa input test for end");
is ($aln->score, 11, 'xmfa alignment score');

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'xmfa');
$status = $strout->write_aln($aln);
is $status, 1,"xmfa output test";

