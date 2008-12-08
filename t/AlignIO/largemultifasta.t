# -*-Perl-*- Test Harness script for Bioperl
# $Id: largemultifasta.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 7);
	
	use_ok('Bio::AlignIO::largemultifasta');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# LARGEMULTIFASTA
$str = Bio::AlignIO->new(
   '-file' => test_input_file('little.largemultifasta'),
                         '-format' => 'largemultifasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Human:/1-81', "fasta input test ";
is ($aln->get_seq_by_pos(1)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    "fasta input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'Rat:',
    "fasta input test for id");

is ($aln->get_seq_by_pos(3)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    "fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(),
                            '-format' => 'largemultifasta');
$status = $strout->write_aln($aln);
is $status, 1,"fasta output test";

