# -*-Perl-*- Test Harness script for Bioperl
# $Id: msf.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 4);
	
	use_ok('Bio::AlignIO::msf');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# MSF
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.msf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', "msf input test";

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			    '-format' => 'msf');
$status = $strout->write_aln($aln);
is $status, 1, "msf output test";

