# -*-Perl-*- Test Harness script for Bioperl
# $Id: mase.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3);
	
	use_ok('Bio::AlignIO::mase');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# MASE
$str = Bio::AlignIO->new(
   '-file' => test_input_file("testaln.mase"),
			   '-format' => 'mase');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/1-318', "mase input test ";

