# -*-Perl-*- Test Harness script for Bioperl
# $Id: xmfa.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 30);
	
	use_ok('Bio::AlignIO::xmfa');
}

my $DEBUG = test_debug(); # foo

my ($str,$aln,$strout,$status);

# XMFA
$str = Bio::AlignIO->new(
		 -file => test_input_file("testaln.xmfa"), 
		 -format => 'xmfa');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');

# test seqs

my @test_data = (
    # 1:1-598 + chrY 
    [ 'chrY/1-598', 1, 598, 1, 'chrY', undef],

    # 2:5000-5534 - chr17 
    [ 'chr17/5534-5000', 5000, 5534, -1, 'chr17', undef],

    # 3:19000-19537 - chr7
    [ 'chr7/19537-19000', 19000, 19537, -1, 'chr7', undef],
);

for my $pos (1..3) {
    my $seq = $aln->get_seq_by_pos($pos);
    my @seq_data = @{shift @test_data};
    is $seq->get_nse, shift @seq_data,  "xmfa input test ";
    is $seq->start, shift @seq_data, "xmfa input test for start";
    is $seq->end, shift @seq_data, "xmfa input test for end";
    is $seq->strand, shift @seq_data,  "xmfa strand test";
    is $seq->display_id, shift @seq_data, "xmfa input test for id";
    is $seq->description, shift @seq_data, "xmfa input test for id";
}

# test aln
is $aln->score, 111, 'xmfa alignment score';

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

