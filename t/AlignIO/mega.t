# -*-Perl-*- Test Harness script for Bioperl
# $Id: mega.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6);
	
	use_ok('Bio::AlignIO::mega');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# MEGA
$str = Bio::AlignIO->new('-format' => 'mega',
  	'-file'   => test_input_file("testaln.mega"));

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'Human/1-141');
is($aln->get_seq_by_pos(2)->get_nse,'Horse/1-144');
$aln->unmatch();
is($aln->get_seq_by_pos(3)->subseq(1,10), 'V-LSAADKGN');

$strout = Bio::AlignIO->new('-format' => 'mega',
	  '-file'   => ">" .test_output_file());

$status = $strout->write_aln($aln);
is $status, 1, "mega output test";


