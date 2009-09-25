# -*-Perl-*- Test Harness script for Bioperl
# $Id: meme.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);
	
	use_ok('Bio::AlignIO::meme');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# MEME
# this file has no Strand column
$str = Bio::AlignIO->new(
		-file => test_input_file('test.meme'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,25;
is $aln->num_sequences,4;
is $aln->get_seq_by_pos(3)->seq(),"CCTTAAAATAAAATCCCCACCACCA";
is $aln->get_seq_by_pos(3)->strand,"1";

# this file has a Strand column
$str = Bio::AlignIO->new(
		-file => test_input_file('test.meme2'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,20;
is $aln->num_sequences,8;
is $aln->get_seq_by_pos(8)->seq(),"CCAGTCTCCCCTGAATACCC";
is $aln->get_seq_by_pos(7)->strand,"-1";
is $aln->get_seq_by_pos(6)->strand,"1";
