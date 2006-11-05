# This is -*-Perl-*- code
# $Id$
use strict;

my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 168;
}
use_ok('Bio::SimpleAlign');
use_ok('Bio::AlignIO');
use_ok('Bio::Root::IO');

END {
    unlink(Bio::Root::IO->catfile("t","data","testout2.pfam"),
	   Bio::Root::IO->catfile("t","data","testout.selex"),
	   Bio::Root::IO->catfile("t","data","testout.pfam"),
	   Bio::Root::IO->catfile("t","data","testout.msf"),
	   Bio::Root::IO->catfile("t","data","testout.fasta"),
	   Bio::Root::IO->catfile("t","data","testout.clustal"),
	   Bio::Root::IO->catfile("t","data","testout.phylip"),
	   Bio::Root::IO->catfile("t","data","testout.nexus"),
	   Bio::Root::IO->catfile("t","data","testout.mega"),
	   Bio::Root::IO->catfile("t","data","testout.po"),
	   Bio::Root::IO->catfile("t","data","testout.largemultifasta")
	  );
}

my ($str,$aln,$strout,$status);

#PSI format  
$str  = new Bio::AlignIO (
    '-file'	=> Bio::Root::IO->catfile("t","data","testaln.psi"),
    '-format'	=> 'psi');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'QUERY/1-798');
is($aln->no_sequences, 56);

# STOCKHOLM (multiple concatenated files, as Pfam flatfile)
$str  = new Bio::AlignIO (
    '-file'	=> Bio::Root::IO->catfile("t","data","testaln.stockholm"),
    '-format'	=> 'stockholm');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');


# PFAM
$str = Bio::AlignIO->new(
	  '-file' => Bio::Root::IO->catfile("t","data","testaln.pfam"));
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');

$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout.pfam"), 
			    '-format' => 'pfam');
$status = $strout->write_aln($aln);
is($status, 1, " failed pfam output test");


# MAF
$str = Bio::AlignIO->new(
	  '-file' => Bio::Root::IO->catfile("t","data","humor.maf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'NM_006987/1-5000', " failed maf input test";
is $aln->get_seq_by_pos(1)->strand, '-';

# MSF
$str = Bio::AlignIO->new(
    '-file' => Bio::Root::IO->catfile("t","data","testaln.msf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed msf input test";

$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout.msf"), 
			    '-format' => 'msf');
$status = $strout->write_aln($aln);
is $status, 1, "  failed msf output test";


# FASTA
$str = Bio::AlignIO->new(
		 -file => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
		 -format => 'fasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', 
  " failed fasta input test ";
is ($aln->get_seq_by_pos(1)->description, 'DESCRIPTION HERE', 
    " failed fasta input test for description");
is ($aln->get_seq_by_pos(11)->display_id, 'AK_YEAST',
    " failed fasta input test for id");

is ($aln->get_seq_by_pos(2)->end, 318,
    " failed fasta input test for end");

is ($aln->get_seq_by_pos(11)->description, 'A COMMENT FOR YEAST', 
    " failed fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout.fasta"), 
			      '-format' => 'fasta');
$status = $strout->write_aln($aln);
is $status, 1,"  failed fasta output test";

my $in = Bio::AlignIO->newFh(
   '-file'  => Bio::Root::IO->catfile("t","data","testaln.fasta"), 
			       '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout2.pfam"), 
				'-format' => 'pfam');
while ( $aln = <$in>) {
    is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431',
     "  failed filehandle input test  ";
    $status = print $out $aln;
    last;
}
is $status, 1, "  failed filehandle output test";


# SELEX
$str = Bio::AlignIO->new(
    '-file' => Bio::Root::IO->catfile("t","data","testaln.selex"),
			   '-format' => 'selex');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', " failed selex format test ";

$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout.selex"), 
			      '-format' => 'selex');
$status = $strout->write_aln($aln);
is $status, 1, "  failed selex output test";


# MASE
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","testaln.mase"),
			   '-format' => 'mase');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/1-318', " failed mase input test ";


# PRODOM
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","testaln.prodom"),
			   '-format' => 'prodom');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', " failed prodom input test ";


# CLUSTAL
$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t","data","testout.clustal"), 
			      '-format' => 'clustalw');
$status = $strout->write_aln($aln);
is $status, 1, "  failed clustalw (.aln) output test";
undef $strout;
$str = Bio::AlignIO->new(
   '-file'=> Bio::Root::IO->catfile("t","data","testout.clustal"), 
			   '-format' => 'clustalw');
$aln = $str->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "  failed clustalw (.aln) input test";
my $io = Bio::AlignIO->new(
   -file => Bio::Root::IO->catfile("t","data","testaln.aln") );
$aln = $io->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string, 'MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEALAQILVPVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEHMDTHYRVLPNIAKLFDDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISNDLKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRVAGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYIDLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI', " failed clustalw consensus_string test";


# BL2SEQ
$str = Bio::AlignIO->new(
   '-file'   => Bio::Root::IO->catfile("t","data","bl2seq.out"),
			 '-format' => 'bl2seq',
			 '-report_type' => 'blastp');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(2)->get_nse, 'ALEU_HORVU/60-360', 
    "failed BLAST bl2seq format test";


# PHYLIP interleaved
$str = Bio::AlignIO->new(
    '-file' => Bio::Root::IO->catfile("t","data","testaln.phylip"),
    '-format' => 'phylip');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new(
    '-file'  => ">".Bio::Root::IO->catfile("t","data","testout.phylip"),
    '-format' => 'phylip');
$status = $strout->write_aln($aln);
is $status, 1, "  failed phylip output test";


# METAFASTA (tests 37-39)
print STDERR "Better Metafasta tests needed\n" if $DEBUG;
$io = Bio::AlignIO->new(-verbose => -1, 
   -file => Bio::Root::IO->catfile("t","data","testaln.metafasta"));
$aln = $io->next_aln;
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string,'CDEFHIJKLMNOPQRSTUVWXYZ', " failed consensus_string on metafasta";
is $aln->symbol_chars,'23'," failed symbol_chars() using metafasta";

# NEXUS
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","testaln.nexus"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapiens/1-45';
$strout = Bio::AlignIO->new('-file'  => ">".
			  Bio::Root::IO->catfile("t", "data", "testout.nexus"),
			    '-format' => 'nexus', );
$status = $strout->write_aln($aln);
is $status, 1, "  failed nexus output test";

$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","Bird_Ovomucoids.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","basic-ladder.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","Kingdoms_DNA.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => Bio::Root::IO->catfile("t","data","char-interleave.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","Primate_mtDNA.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","char-matrix-spaces.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","SPAN_Family4nl.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","SPAN_Family7n.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","long-names.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","SPAN_Family8a.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","multiline-intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => Bio::Root::IO->catfile("t","data","Treebase-chlamy-dna.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => Bio::Root::IO->catfile("t","data","quoted-strings1.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","UnaSmithHIV-both.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => Bio::Root::IO->catfile("t","data","quoted-strings2.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","barns-combined.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","radical-whitespace.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t","data","basic-bush.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => Bio::Root::IO->catfile("t","data","radical-whitespace_02.nex"),
			  '-format' => 'nexus');


# EMBOSS water
$str = new Bio::AlignIO('-format' => 'emboss',
		 '-file'   => Bio::Root::IO->catfile("t", "data", 'cysprot.water'));
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
$str = new Bio::AlignIO('-format' => 'emboss',
	  '-file'   => Bio::Root::IO->catfile("t", "data", 'cysprot.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'499.50');
is($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/1-345');
is($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-333');


# EMBOSS water 2.2.x
$str = new Bio::AlignIO('-format' => 'emboss',
	 '-file'   => Bio::Root::IO->catfile("t", "data", 'cys1_dicdi.water'));
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
$str = new Bio::AlignIO(-verbose => $DEBUG,
	  '-format' => 'emboss',
   	'-file'   => Bio::Root::IO->catfile("t", "data", 'sparsealn.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'18.0');
is(sprintf("%.1f",$aln->overall_percentage_identity), 2.1);
is(sprintf("%.1f",$aln->average_percentage_identity), 38.5);
is($aln->get_seq_by_pos(1)->length, 238);
is($aln->length,238);
is($aln->get_seq_by_pos(1)->get_nse,'KV1K_HUMAN/1-108');
is($aln->get_seq_by_pos(2)->get_nse,'IF1Y_HUMAN/1-143');
is($aln->get_seq_by_pos(1)->seq(), 'DIQMTQSPSTLSVSVGDRVTITCEASQTVLSYLNWYQQKPGKAPKLLIYAASSLETGVPSRFSGQGSGTBFTFTISSVZPZBFATYYCQZYLDLPRTFGQGTKVDLKR'.'-'x130);
is($aln->get_seq_by_pos(2)->seq(), ('-'x94).'PKNKGKGGK-NRRRGKNENESEKRELVFKEDGQEYAQVIKMLGNGRLEALCFDGVKRLCHIRGKLRKKVWINTSDIILVGLRDYQDNKADVILKYNADEARSLKAYGGLPEHAKINETDTFGPGDDDEIQFDDIGDDDEDIDDI');
is($aln->is_flush, 1);


# MEGA
$str = new Bio::AlignIO('-format' => 'mega',
  	'-file'   => Bio::Root::IO->catfile("t","data","hemoglobinA.meg"));

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'Human/1-141');
is($aln->get_seq_by_pos(2)->get_nse,'Horse/1-144');
$aln->unmatch();
is($aln->get_seq_by_pos(3)->subseq(1,10), 'V-LSAADKGN');

$strout = new Bio::AlignIO('-format' => 'mega',
	  '-file'   => ">" .Bio::Root::IO->catfile("t","data","testout.mega"));

$status = $strout->write_aln($aln);
is $status, 1, "  failed mega output test";


# EMBOSS needle
$str = new Bio::AlignIO('-format' => 'emboss',
	  '-file'   => Bio::Root::IO->catfile('t','data','gf-s71.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'MEDVTLFQFTWRKPI-RLQGEIVYKTSETQTIETNKKDVECVANFQENKEVQTDS-VDNGVGENVKKDITISKEVLNLLYDFVRDDSKVNYDRLLEFHKFDKVALETVQKYHVETRNENIILMISSSSRKTLILFGGISHETFCSHQARALLCSSSTSFSIPLPVCAISAVFYSSTQFILGDVSGNISMCSKDKIIFEKKITDGAVTCLEMCRHGLLSGSDDGNIILWQIGTSGLEKLGGTKLTVSDLSRKIRRSSTSNKPVAIVSMQVYVWPSGEEACVATETGGLYLLTLPTLDYKPLSHQTATSINKILFENQFVAVIYHTSNAAVFNSEGLVDEIPFVATLAVR----------PKLVLF--YTSVCVQDITLNCTSPFREFNNEYNPVIKFSKIRFSADLSVING-FRTSSPNSNN-----------------------------------------------');
is($aln->get_seq_by_pos(1)->seq(), 'MEDVTLHHFRWRKPVENKNGEIVYKTSETQTAEISRKDVECVANFQKSQESQTDDFMQNGVGDGIKKEIRISKEVLGHIYDFLRDDSKVNYDRLLEFHKFDKVSLETVQKYHVETRNENIILMISNSSRKTLILFGGLSHETFCSHQARAVLCSSSTTSSLPLPVCAISAVFYSSTQFLLGDISGNISMWTKEKMIFENKVTDGSVTSLELCRYGLLSGSDDGNVILWKVEESKIEKIEGIKLTVSDLSRKIRRSSTSNKPVAIVSMQV----SGDEVCVATETGGLYLLTLPTLESKPLT-QSATSIFKILYEHPYIAVVYHTSNSAIFNSEGLVDEIPFVATLAVRCGAYFIFSNQSRLIIWSMNTRSTVIDENLNCHS-ICSLSND--------------TLQVLDGDFNLNSQSENSATSESENLRISDLQNLRMLKLQNLRTSEFQNFRTSESQYFKKDNGEL');
is($aln->is_flush(), 1);
is($aln->get_seq_by_pos(1)->get_nse,'gf.s71.44/1-448');
is($aln->get_seq_by_pos(2)->get_nse,'Y50C1A.2/1-406');


# PHYLIP sequential/non-interleaved
$strout = Bio::AlignIO->new('-file'  =>
			    Bio::Root::IO->catfile("t","data",
						   "noninterleaved.phy"),
			    '-format' => 'phylip');
$aln = $strout->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAAAGGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGAATTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGGTTTATCAAAGTAAGACAGTATGATCAGATACCCATAGAGATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCCACACCTGTCAATATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' );


# LARGEMULTIFASTA
$str = Bio::AlignIO->new(
   '-file' => Bio::Root::IO->catfile("t", "data","little.largemultifasta"),
                         '-format' => 'largemultifasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Human:/1-81', " failed fasta input test ";
is ($aln->get_seq_by_pos(1)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    " failed fasta input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'Rat:',
    " failed fasta input test for id");

is ($aln->get_seq_by_pos(3)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    " failed fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".Bio::Root::IO->catfile("t", "data",
                                       "testout.largemultifasta"),
                            '-format' => 'largemultifasta');
$status = $strout->write_aln($aln);
is $status, 1,"  failed fasta output test";


# POA
# just skip on perl 5.6.0 and earlier as it causes a crash on 
# default perl with OS X 10.2
# fink perl 5.6.0 does not seem to have the problem
# can't figure out what it is so just skip for now
if( $^O ne 'darwin' || $] > 5.006 ) {
	$str = new Bio::AlignIO(
			  -file   => Bio::Root::IO->catfile("t", "data", "testaln.po"),
			  -format => 'po',
			  );
	isa_ok($str, 'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	is $aln->no_sequences, 6;

# output is? i.e. does conversion from clustalw to po give the same alignment?
	$str = new Bio::AlignIO(
		  '-file'   => Bio::Root::IO->catfile("t", "data", "testaln.aln"),
		  '-format' => 'clustalw');
	isa_ok($str,'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	$strout = Bio::AlignIO->new(
		 '-file'   => ">" . Bio::Root::IO->catfile("t", "data", "testout.po"),
		 '-format' => 'po');
	$status = $strout->write_aln($aln);
	is $status, 1, " failed po output test";

	$str = new Bio::AlignIO(
		 '-file'   => Bio::Root::IO->catfile("t", "data", "testaln.po"),
		 '-format' => 'po');
	isa_ok($str,'Bio::AlignIO');
	my $aln2 = $str->next_aln();
	isa_ok($aln2,'Bio::Align::AlignI');
	is $aln2->no_sequences, $aln->no_sequences;
	is $aln2->length, $aln->length;
} else {
	for ( 1..14 ) {
		skip(1,"skipping due to bug in perl 5.6.0 that comes with OS X 10.2");
	}
}


# MEME
# this file has no Strand column
$str = new Bio::AlignIO(
		-file => Bio::Root::IO->catfile("t", "data", "test.meme"),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,25;
is $aln->no_sequences,4;
is $aln->get_seq_by_pos(3)->seq(),"CCTTAAAATAAAATCCCCACCACCA";
is $aln->get_seq_by_pos(3)->strand,"1";

# this file has a Strand column
$str = new Bio::AlignIO(
		-file => Bio::Root::IO->catfile("t", "data", "test.meme2"),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,20;
is $aln->no_sequences,8;
is $aln->get_seq_by_pos(8)->seq(),"CCAGTCTCCCCTGAATACCC";
is $aln->get_seq_by_pos(7)->strand,"-1";
is $aln->get_seq_by_pos(6)->strand,"1";
