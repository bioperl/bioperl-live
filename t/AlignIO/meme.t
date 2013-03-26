# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 20);
	
	use_ok('Bio::AlignIO::meme');
}

my $DEBUG = test_debug();

# MEME
# this file has no Strand column, and it's version 3.0
my $str = Bio::AlignIO->new(
		-file => test_input_file('test-3.0-1.meme'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
my $aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,25;
is $aln->num_sequences,4;
is $aln->get_seq_by_pos(3)->seq(),"CCTTAAAATAAAATCCCCACCACCA";
is $aln->get_seq_by_pos(3)->strand,"1";

# this file has a Strand column, also version 3.0
$str = Bio::AlignIO->new(
		-file => test_input_file('test-3.0-2.meme'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,20;
is $aln->num_sequences,8;
is $aln->get_seq_by_pos(8)->seq(),"CCAGTCTCCCCTGAATACCC";
is $aln->get_seq_by_pos(7)->strand,"-1";
is $aln->get_seq_by_pos(6)->strand,"1";

# version 4.9
$str = Bio::AlignIO->new(
		-file => test_input_file('test-4.9.meme'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,21;
is $aln->num_sequences,47;
is $aln->get_seq_by_pos(3)->seq(),"AGAGAAACAAGAGGCCTCTTT";
is $aln->get_seq_by_pos(3)->strand,"1";
