# This is -*-Perl-*- code
# $Id$
use strict;

my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN { 
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 80;
}

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::IO;

ok (1);

END {

    unlink(Bio::Root::IO->catfile("t","data","testout2.pfam"),
 	   Bio::Root::IO->catfile("t","data","testout.selex"),
 	   Bio::Root::IO->catfile("t","data","testout.pfam"),
 	   Bio::Root::IO->catfile("t","data","testout.msf"),
 	   Bio::Root::IO->catfile("t","data","testout.fasta"), 
 	   Bio::Root::IO->catfile("t","data","testout.clustal"),
 	   Bio::Root::IO->catfile("t","data","testout.phylip"),
 	   Bio::Root::IO->catfile("t","data","testout.nexus"),
 	   Bio::Root::IO->catfile("t","data","testout.mega"));
    unlink(
           Bio::Root::IO->catfile("t","data","littleout.largemultifasta")
          );
}

my ($str,$aln,$strout,$status);

# STOCKHOLM (multiple concatenated files, as Pfam flatfile)
ok $str  = new Bio::AlignIO (
    '-file'	=> Bio::Root::IO->catfile("t","data","testaln.stockholm"),
    '-format'	=> 'stockholm');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246';
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246';
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246';

# PFAM
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed pfam input test";

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.pfam"), 
			      '-format' => 'pfam');
$status = $strout->write_aln($aln);
ok $status, 1, " failed pfam output test";

# MAF
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","humor.maf"));
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 
    'NM_006987/0-5000', " failed maf input test";

# MSF
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.msf"));
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed msf input test";


$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.msf"), 
			    '-format' => 'msf');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed msf output test";

# FASTA
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
			   '-format' => 'fasta');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', " failed fasta input test ";
ok ($aln->get_seq_by_pos(1)->description, 'DESCRIPTION HERE', 
    " failed fasta input test for description");
ok ($aln->get_seq_by_pos(11)->display_id, 'AK_YEAST',
    " failed fasta input test for id");

ok ($aln->get_seq_by_pos(11)->description, 'A COMMENT FOR YEAST', 
    " failed fasta input test for description");

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.fasta"), 
			      '-format' => 'fasta');
$status = $strout->write_aln($aln);
ok $status, 1,"  failed fasta output test";

# SELEX
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.selex"),
			   '-format' => 'selex');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', " failed selex format test ";

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.selex"), 
			      '-format' => 'selex');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed selex output test";

# MASE
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.mase"),
			   '-format' => 'mase');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/1-318', " failed mase input test ";

# PRODOM
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.prodom"),
			   '-format' => 'prodom');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', " failed prodom input test ";

# CLUSTAL
$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t","data","testout.clustal"), 
			      '-format' => 'clustalw');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed clustalw (.aln) output test";
undef $strout;
$str = Bio::AlignIO->new('-file'=> Bio::Root::IO->catfile("t","data","testout.clustal"), 
			   '-format' => 'clustalw');
$aln = $str->next_aln($aln);
ok $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "  failed clustalw (.aln) input test";
my $io = Bio::AlignIO->new(-file => Bio::Root::IO->catfile("t","data","testaln.aln") );
$aln = $io->next_aln();
ok $aln->consensus_string, 'MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEALAQILVPVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEHMDTHYRVLPNIAKLFDDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI', " failed clustalw consensus_string test";

# FASTA
my $in = Bio::AlignIO->newFh('-file'  => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
			       '-format' => 'fasta');
my $out = Bio::AlignIO->newFh('-file'  => ">".Bio::Root::IO->catfile("t","data","testout2.pfam"), 
				'-format' => 'pfam');
while ( $aln = <$in>) {
    ok $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', "  failed filehandle input test  ";
    $status = print $out $aln;
    last;
}
ok $status, 1, "  failed filehandle output test";

# BL2SEQ
$str = Bio::AlignIO->new('-file'   => Bio::Root::IO->catfile("t","data","bl2seq.out"),
			 '-format' => 'bl2seq',
			 '-report_type' => 'blastp');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(2)->get_nse, 'ALEU_HORVU/60-360', 
    "failed BLAST bl2seq format test";

# PHYLIP
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.phylip"),
			  '-format' => 'phylip');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new('-file'  => ">".Bio::Root::IO->catfile("t","data","testout.phylip"),
			      '-format' => 'phylip');
$status = $strout->write_aln($aln);
ok $status, 1, "  failed phylip output test";

# METAFASTA

$io = Bio::AlignIO->new(-file => Bio::Root::IO->catfile("t","data","testaln.metafasta"));
$aln = $io->next_aln;
ok $aln->consensus_string,'CDEFHIJKLMNOPQRSTUVWXYZhydrophobicIOIOIJOIIOOIOOOOUIIXstructuralABAEEIEIJEIIEOAEEAAUIAX', " failed consensus_string on metafasta";
ok $aln->percentage_identity,'100', " failed percentage_identity using metafasta";
ok $aln->symbol_chars,'39'," failed symbol_chars() using metafasta";

# NEXUS
$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t","data","testaln.nexus"),
			  '-format' => 'nexus');
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapiens/1-45';
$strout = Bio::AlignIO->new('-file'  => ">".
			  Bio::Root::IO->catfile("t", "data", "testout.nexus"),
			    '-format' => 'nexus', );
$status = $strout->write_aln($aln);
ok $status, 1, "  failed nexus output test";

# EMBOSS water

$str = new Bio::AlignIO('-format' => 'emboss',
			'-file'   => Bio::Root::IO->catfile("t", "data", 'cysprot.water'));
$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/3-342');
ok($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-331');
ok(sprintf("%.1f",$aln->overall_percentage_identity),33.8);
ok(sprintf("%.1f",$aln->average_percentage_identity),40.1);

ok($aln->get_seq_by_pos(1)->start, 3);
ok($aln->length,364);

# EMBOSS needle

$str = new Bio::AlignIO('-format' => 'emboss',
			'-file'   => Bio::Root::IO->catfile("t", "data", 'cysprot.needle'));
$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/1-345');
ok($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-333');


# EMBOSS water 2.2.x

$str = new Bio::AlignIO('-format' => 'emboss',
			'-file'   => Bio::Root::IO->catfile("t", "data", 'cys1_dicdi.water'));
$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/1-343');
ok($aln->get_seq_by_pos(2)->get_nse,'CYS1_DICDI-1/1-343');

$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/29-343');
ok($aln->get_seq_by_pos(2)->get_nse,'ALEU_HORVU/61-360');

# EMBOSS water 2.2.x sparse needle

$str = new Bio::AlignIO(-verbose => $DEBUG,
			'-format' => 'emboss',
			'-file'   => Bio::Root::IO->catfile("t", "data", 'sparsealn.needle'));
$aln = $str->next_aln();
ok($aln);
ok(sprintf("%.1f",$aln->overall_percentage_identity), 2.1);
ok(sprintf("%.1f",$aln->average_percentage_identity), 38.5);
ok($aln->get_seq_by_pos(1)->length, 238);
ok($aln->length,238);
ok($aln->get_seq_by_pos(1)->get_nse,'KV1K_HUMAN/1-108');
ok($aln->get_seq_by_pos(2)->get_nse,'IF1Y_HUMAN/1-143');
ok($aln->get_seq_by_pos(1)->seq(), 'DIQMTQSPSTLSVSVGDRVTITCEASQTVLSYLNWYQQKPGKAPKLLIYAASSLETGVPSRFSGQGSGTBFTFTISSVZPZBFATYYCQZYLDLPRTFGQGTKVDLKR'.'-'x130);
ok($aln->get_seq_by_pos(2)->seq(), ('-'x94).'PKNKGKGGK-NRRRGKNENESEKRELVFKEDGQEYAQVIKMLGNGRLEALCFDGVKRLCHIRGKLRKKVWINTSDIILVGLRDYQDNKADVILKYNADEARSLKAYGGLPEHAKINETDTFGPGDDDEIQFDDIGDDDEDIDDI');
ok($aln->is_flush);
# MEGA

$str = new Bio::AlignIO('-format' => 'mega',
			'-file'   => Bio::Root::IO->catfile("t","data","hemoglobinA.meg"));

$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(1)->get_nse,'Human/1-141');
ok($aln->get_seq_by_pos(2)->get_nse,'Horse/1-144');
$aln->unmatch();
ok($aln->get_seq_by_pos(3)->subseq(1,10), 'V-LSAADKGN');

$strout = new Bio::AlignIO('-format' => 'mega',
			   '-file'   => ">" .Bio::Root::IO->catfile("t","data","testout.mega"));

$status = $strout->write_aln($aln);
ok $status, 1, "  failed mega output test";


# 

$str = new Bio::AlignIO('-format' => 'emboss',
			'-file'   => Bio::Root::IO->catfile('t','data','gf-s71.needle'));
$aln = $str->next_aln();
ok($aln);
ok($aln->get_seq_by_pos(2)->seq(), 'MEDVTLFQFTWRKPI-RLQGEIVYKTSETQTIETNKKDVECVANFQENKEVQTDS-VDNGVGENVKKDITISKEVLNLLYDFVRDDSKVNYDRLLEFHKFDKVALETVQKYHVETRNENIILMISSSSRKTLILFGGISHETFCSHQARALLCSSSTSFSIPLPVCAISAVFYSSTQFILGDVSGNISMCSKDKIIFEKKITDGAVTCLEMCRHGLLSGSDDGNIILWQIGTSGLEKLGGTKLTVSDLSRKIRRSSTSNKPVAIVSMQVYVWPSGEEACVATETGGLYLLTLPTLDYKPLSHQTATSINKILFENQFVAVIYHTSNAAVFNSEGLVDEIPFVATLAVR----------PKLVLF--YTSVCVQDITLNCTSPFREFNNEYNPVIKFSKIRFSADLSVING-FRTSSPNSNN-----------------------------------------------');
ok($aln->get_seq_by_pos(1)->seq(), 'MEDVTLHHFRWRKPVENKNGEIVYKTSETQTAEISRKDVECVANFQKSQESQTDDFMQNGVGDGIKKEIRISKEVLGHIYDFLRDDSKVNYDRLLEFHKFDKVSLETVQKYHVETRNENIILMISNSSRKTLILFGGLSHETFCSHQARAVLCSSSTTSSLPLPVCAISAVFYSSTQFLLGDISGNISMWTKEKMIFENKVTDGSVTSLELCRYGLLSGSDDGNVILWKVEESKIEKIEGIKLTVSDLSRKIRRSSTSNKPVAIVSMQV----SGDEVCVATETGGLYLLTLPTLESKPLT-QSATSIFKILYEHPYIAVVYHTSNSAIFNSEGLVDEIPFVATLAVRCGAYFIFSNQSRLIIWSMNTRSTVIDENLNCHS-ICSLSND--------------TLQVLDGDFNLNSQSENSATSESENLRISDLQNLRMLKLQNLRTSEFQNFRTSESQYFKKDNGEL');
ok($aln->is_flush());
ok($aln->get_seq_by_pos(1)->get_nse,'gf.s71.44/1-448');
ok($aln->get_seq_by_pos(2)->get_nse,'Y50C1A.2/1-406');


# Let's parse Phylip non-interleaved
$strout = Bio::AlignIO->new('-file'  =>
			    Bio::Root::IO->catfile("t","data",
						   "noninterleaved.phy"),
			    '-format' => 'phylip');
$aln = $strout->next_aln($aln);
ok($aln);
ok($aln->get_seq_by_pos(2)->seq(), 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGGTTTATCAAAGTAAGACAGTATGATCAGATACCCATAGAGATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCCACACCTGTCAATATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' );

#LARGEMULTIFASTA

$str = Bio::AlignIO->new('-file' => Bio::Root::IO->catfile("t", "data","little.largemultifasta"),
                         '-format' => 'largemultifasta');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, 'Human:/1-81', " failed fasta input test ";
ok ($aln->get_seq_by_pos(1)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    " failed fasta input test for description");
ok ($aln->get_seq_by_pos(3)->display_id, 'Rat:',
    " failed fasta input test for id");

ok ($aln->get_seq_by_pos(3)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    " failed fasta input test for description");

$strout = Bio::AlignIO->new('-file' => ">".Bio::Root::IO->catfile("t", "data",
                                                                  "littleout.largemultifasta"),
                            '-format' => 'largemultifasta');
$status = $strout->write_aln($aln);
ok $status, 1,"  failed fasta output test";
