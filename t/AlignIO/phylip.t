# -*-Perl-*- Test Harness script for Bioperl
# $Id: phylip.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 17);

	use_ok('Bio::AlignIO::phylip');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PHYLIP sequential/non-interleaved
$strout = Bio::AlignIO->new('-file'  => test_input_file('noninterleaved.phy'), '-interleaved' => 0,
			    '-format' => 'phylip');
$aln = $strout->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAA'.
   'AGGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGAATT'.
   'TGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGGTTTATCAAAGTAAGACAGTATGATCAGA'.
   'TACCCATAGAGATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCCACACCTGTCAATATAATTG'.
   'GAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' );

# PHYLIP interleaved with long Ids

$str = Bio::AlignIO->new(
    '-file' => test_input_file("protpars_longid.phy"),
    '-format' => 'phylip',
    'longid' => 1);

$aln = $str->next_aln();
#isa_ok($str,'Bio::AlignIO');
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'S I N F R U  P 0 0 1 /1-84';
is $aln->get_seq_by_pos(2)->get_nse, 'SINFRUP002/1-84';

# PHYLIP interleaved, multiple segments
$str = Bio::AlignIO->new(
    '-file' => test_input_file("protpars.phy"),
    '-format' => 'phylip');

$aln = $str->next_aln();
#isa_ok($str,'Bio::AlignIO');
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'SINFRUP001/1-4940';
# is $aln->get_seq_by_pos(2)->get_nse, 'SINFRUP002/1-84';


# PHYLIP interleaved

$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.phylip"),
    '-format' => 'phylip');
$aln = $str->next_aln();
#isa_ok($str,'Bio::AlignIO');
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new(
    '-file'  => ">".test_output_file(),
    '-format' => 'phylip');
$status = $strout->write_aln($aln);
is $status, 1, "phylip output test";

# check the LocatableSeq start/end/strand etc
my $ls = $aln->get_seq_by_pos(2);
is($ls->display_id, 'Pan_panisc');
is($ls->start, 1);
is($ls->end,47);

# bug 2984
TODO: {
    local $TODO = 'problems with default strand, length?';
    # shouldn't this be 0?
    is($ls->strand,0);
    is($ls->length,47);
}

# check to see that newlines between header and sequences are parsed correctly
$str = Bio::AlignIO->new('-file' => test_input_file("codeml45b.mlc"), '-format' => 'phylip', '-longid' => 1);
$aln = $str->next_aln();
$ls = $aln->get_seq_by_pos(9);
ok($ls->display_id eq "Pop_trich_ch", "newline between header and sequences is parsed correctly");
