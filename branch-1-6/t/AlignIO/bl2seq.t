# -*-Perl-*- Test Harness script for Bioperl
# $Id: bl2seq.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3);
	
	use_ok('Bio::AlignIO::bl2seq');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# BL2SEQ
$str = Bio::AlignIO->new(
   '-file'   => test_input_file("bl2seq.out"),
			 '-format' => 'bl2seq',
			 '-report_type' => 'blastp');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(2)->get_nse, 'ALEU_HORVU/60-360', 
    "BLAST bl2seq format test";