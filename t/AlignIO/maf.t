# -*-Perl-*- Test Harness script for Bioperl
# $Id: maf.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 11);
	
	use_ok('Bio::AlignIO::maf');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# MAF
$str = Bio::AlignIO->new(
	  '-file' => test_input_file("humor.maf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'NM_006987/5000-1', "maf input test";
is $aln->get_seq_by_pos(1)->strand, -1;

# MAF - bug 2453
$str = Bio::AlignIO->new(
	  '-file' => test_input_file("bug2453.maf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'hg16.chr7/27578829-27578866', "maf input test";
is $aln->get_seq_by_pos(1)->strand, 1;
$aln = $str->next_aln();
is $aln->get_seq_by_pos(1)->get_nse, 'hg16.chr7/27699740-27699745', "maf input test";
is $aln->get_seq_by_pos(1)->strand, 1;
$aln = $str->next_aln();
is $aln->get_seq_by_pos(1)->get_nse, 'hg16.chr7/27707222-27707234', "maf input test";
is $aln->get_seq_by_pos(1)->strand, 1;

