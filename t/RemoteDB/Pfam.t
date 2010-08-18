# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

#test Bio::DB::Align::Pfam
BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 35,
			   -requires_modules => [qw(HTTP::Request
									    LWP::UserAgent
										Bio::AlignIO
										Bio::DB::GenericWebAgent
										Bio::DB::Align)],
			   -requires_networking => 1);
	
	use_ok('Bio::Root::IO');
	use_ok('Bio::Root::Root');
}


#test calling from Bio::DB::Align
my ($dbobj, $aln);
ok $dbobj=Bio::DB::Align->new(-db=>"Pfam"), 'Bio::DB::Align::Pfam';
lives_ok{$aln=$dbobj->get_Aln_by_id("Piwi")};

is $aln->id,"Piwi";
is $aln->accession,"PF02171";
is $aln->length,"412";
is $aln->num_sequences,"21";
is $aln->source,"Pfam";

my @names;
foreach my $seq ($aln->next_Seq) {
	push @names,$seq->display_id;
}
is join(";",@names),"YQ53_CAEEL;Q21691_CAEEL;O48771_ARATH;Q9ZVD5_ARATH;TAG76_CAEEL;O16720_CAEEL;PINH_ARATH;AGO1_SCHPO;O76922_DROME;PIWI_DROME;Q17567_CAEEL;PIWL1_HUMAN;PIWI_ARCFU;Y1321_METJA;O67434_AQUAE;Q21495_CAEEL;O16386_CAEEL;O02095_CAEEL;Q19645_CAEEL;Q23415_CAEEL;O62275_CAEEL";

my $acc=$dbobj->id2acc("Piwi");
is $acc,"PF02171";
my $id=$dbobj->acc2id("PF02171");
is $id,"Piwi";

#Test parameter based calling
my ($dbobj2,$aln2);
ok $dbobj2=Bio::DB::Align::Pfam->new(), 'Bio::DB::Align::Pfam';
lives_ok{$aln2=$dbobj2->get_Aln_by_acc(-accession=>"PF02171",-alignment=>"full",-format=>"stockholm",-order=>"a",-case=>"u",-gap=>"dots")};
is $aln->id,"Piwi";
is $aln->accession,"PF02171";
is $aln2->length,"1030";
is $aln2->num_sequences,"756";
is $aln2->gap_char,".";
is $aln2->source,"stockholm";

my $seq=$aln2->get_Seq_by_pos(1);
is $seq->seq,"............................................................................................L.IVTFINE..KDKDR.IYGQMKKY.CF...QE..H..GISHQNILS.K..F.LKSKNP.................SSVASKIAQ.......QMSMK..........LGNP......LWA.IPKP.NG........I......................S.....D.K..T.M..................................VIGIDIYHKLLTN.................RRSCMGFVAYLE...S...E.C.L...N....T.F.....ARP.I...IM.....K..........E.......G....Q...EM......CHE....VGR...................................................I.TVEAISA.Y...........FE......R.....N...........G........K.........K..YLPD...TIIVFR...............DGVGNAQIEALKQTEILQM.KNA.IR........................SIN...KNYN......PQFAVI.MINK..KIN..DR.......F.FMVNGGGGQN..............................QQKQTLSNP....P.S......GSVI.........AD.KI.T.....SS.......N......F...D..Y.FIT......A......Q..Y.V..T......Q.....GTC...TPTHYRVL..ENDT..................NW.S.EE...L......F....WQF...............TYY.Q.CFNY.QNW.T....GAVRV..................PSCVQYAHKLAYLIGDTY.Q..................................................................";
is $seq->display_id,"A0CB11_PARTE";

#Test parametr based calling from get_Aln_by_id
my $aln3;
lives_ok{$aln3=$dbobj->get_Aln_by_id(-id=>"Kazal_1",-alignment=>"seed",-format=>"msf",-gap=>"dashes")};
is $aln3->id,"Kazal_1";
is $aln3->accession,"PF00050";
is $aln3->length,"66";
is $aln3->num_sequences,"51"; 
is $aln3->gap_char,"-";
my $seq2=$aln3->get_Seq_by_pos(1);
is $seq2->seq,"CKVYTEA----------CT---REYNPICDSAAKTYSNECTF----CNEKM-NNDADIHFNHFGEC";
is $seq2->get_nse,"IAC1_BOVIN/14-61";

#Test selex format
my $aln4=$dbobj->get_Aln_by_id("Kazal_2","seed","selex");
is $aln4->id,"Kazal_2";
is $aln4->accession,"PF07648";
is $aln4->num_sequences,"74";

my $seq3=$aln4->get_Seq_by_pos(1);
is $seq3->seq,"CKPHKR--P--VC------GSNGK---TYLN--HCELHRDACLT-----------GSK--IQVDYDGH-C";
is $seq3->get_nse,"FSTL1_MOUSE/56-96";
