# -*-Perl-*- Test Harness script for Bioperl
# $Id: emboss.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 37);
	
	use_ok('Bio::AlignIO::emboss');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# EMBOSS water
$str = Bio::AlignIO->new('-format' => 'emboss',
		 '-file'   => test_input_file('cysprot.water'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'501.50');
is($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/3-342');
is($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-331');
is(sprintf("%.1f",$aln->overall_percentage_identity),33.8);
is(sprintf("%.1f",$aln->average_percentage_identity),40.1);

is($aln->get_seq_by_pos(1)->start, 3);
is($aln->length,364);


# EMBOSS needle
$str = Bio::AlignIO->new('-format' => 'emboss',
	  '-file'   => test_input_file('cysprot.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'499.50');
is($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/1-345');
is($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-333');


# EMBOSS water 2.2.x
$str = Bio::AlignIO->new('-format' => 'emboss',
	 '-file'   => test_input_file('cys1_dicdi.water'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/1-343');
is($aln->get_seq_by_pos(2)->get_nse,'CYS1_DICDI-1/1-343');
is($aln->score,'1841.0');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/29-343');
is($aln->get_seq_by_pos(2)->get_nse,'ALEU_HORVU/61-360');


# EMBOSS water 2.2.x sparse needle
$str = Bio::AlignIO->new(-verbose => $DEBUG,
	  '-format' => 'emboss',
   	'-file'   => test_input_file('sparsealn.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'18.0');
is(sprintf("%.1f",$aln->overall_percentage_identity), 2.1);
is(sprintf("%.1f",$aln->average_percentage_identity), 38.5);
is($aln->get_seq_by_pos(1)->length, 238);
is($aln->length,238);
is($aln->get_seq_by_pos(1)->get_nse,'KV1K_HUMAN/1-108');
is($aln->get_seq_by_pos(2)->get_nse,'IF1Y_HUMAN/1-143');
is($aln->get_seq_by_pos(1)->seq(), 'DIQMTQSPSTLSVSVGDRVTITCEASQTVLSYLNWYQQK'.
   'PGKAPKLLIYAASSLETGVPSRFSGQGSGTBFTFTISSVZPZBFATYYCQZYLDLPRTFGQGTKVDLKR'.
   '-'x130);
is($aln->get_seq_by_pos(2)->seq(), ('-'x94).'PKNKGKGGK-NRRRGKNENESEKRELVFKE'.
   'DGQEYAQVIKMLGNGRLEALCFDGVKRLCHIRGKLRKKVWINTSDIILVGLRDYQDNKADVILKYNADEAR'.
   'SLKAYGGLPEHAKINETDTFGPGDDDEIQFDDIGDDDEDIDDI');
is($aln->is_flush, 1);

# EMBOSS needle
$str = Bio::AlignIO->new('-format' => 'emboss',
	  '-file'   => test_input_file('gf-s71.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'MEDVTLFQFTWRKPI-RLQGEIVYKTSETQTIETNKKDVECVANFQENKEVQTDS-VDNGVGENVKKDITISKEVLNLLYDFVRDDSKVNYDRLLEFHKFDKVALETVQKYHVETRNENIILMISSSSRKTLILFGGISHETFCSHQARALLCSSSTSFSIPLPVCAISAVFYSSTQFILGDVSGNISMCSKDKIIFEKKITDGAVTCLEMCRHGLLSGSDDGNIILWQIGTSGLEKLGGTKLTVSDLSRKIRRSSTSNKPVAIVSMQVYVWPSGEEACVATETGGLYLLTLPTLDYKPLSHQTATSINKILFENQFVAVIYHTSNAAVFNSEGLVDEIPFVATLAVR----------PKLVLF--YTSVCVQDITLNCTSPFREFNNEYNPVIKFSKIRFSADLSVING-FRTSSPNSNN-----------------------------------------------');
is($aln->get_seq_by_pos(1)->seq(), 'MEDVTLHHFRWRKPVENKNGEIVYKTSETQTAEISRKDVECVANFQKSQESQTDDFMQNGVGDGIKKEIRISKEVLGHIYDFLRDDSKVNYDRLLEFHKFDKVSLETVQKYHVETRNENIILMISNSSRKTLILFGGLSHETFCSHQARAVLCSSSTTSSLPLPVCAISAVFYSSTQFLLGDISGNISMWTKEKMIFENKVTDGSVTSLELCRYGLLSGSDDGNVILWKVEESKIEKIEGIKLTVSDLSRKIRRSSTSNKPVAIVSMQV----SGDEVCVATETGGLYLLTLPTLESKPLT-QSATSIFKILYEHPYIAVVYHTSNSAIFNSEGLVDEIPFVATLAVRCGAYFIFSNQSRLIIWSMNTRSTVIDENLNCHS-ICSLSND--------------TLQVLDGDFNLNSQSENSATSESENLRISDLQNLRMLKLQNLRTSEFQNFRTSESQYFKKDNGEL');
is($aln->is_flush(), 1);
is($aln->get_seq_by_pos(1)->get_nse,'gf.s71.44/1-448');
is($aln->get_seq_by_pos(2)->get_nse,'Y50C1A.2/1-406');

