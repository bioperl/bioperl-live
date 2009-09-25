# -*-Perl-*- Test Harness script for Bioperl
# $Id: prodom.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3);
	
	use_ok('Bio::AlignIO::prodom');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PRODOM
$str = Bio::AlignIO->new(
   '-file' => test_input_file("testaln.prodom"),
			   '-format' => 'prodom');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "prodom input test ";

