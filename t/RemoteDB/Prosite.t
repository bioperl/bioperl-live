# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

#test Bio::DB::Align::Pfam
BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 16,
			   -requires_modules => [qw(HTTP::Request
									    LWP::UserAgent
										Bio::AlignIO
										Bio::DB::GenericWebAgent
										Bio::DB::Align)],
			   -requires_networking => 1);
	
	use_ok('Bio::Root::IO');
	use_ok('Bio::Root::Root');
}

#test for pfam format is not included here, because this format is not well supported by Bio::AlignIO::Pfam

#test calling from Bio::DB::Align
my ($dbobj, $aln);
SKIP: {
	ok $dbobj=Bio::DB::Align->new(-db=>"Prosite"), 'Bio::DB::Align::Prosite';
	eval{$aln=$dbobj->get_Aln_by_acc("PS51092")};
	skip("Bio::DB::Align::Prosite HTTP error:$@", 5) if $@;
	is $aln->accession,"PS51092";
	is $aln->length,"49";
	is $aln->num_sequences,"103";
	is $aln->source,"Prosite";
	
	my @names;
	foreach my $seq ($aln->next_Seq) {
		push @names,$seq->display_id;
	}
	is join(";",@names),"BSPH1_HUMAN;BSPH1_HUMAN;BSPH1_MOUSE;BSPH1_MOUSE;ESPB1_CANFA;ESPB1_CANFA;ESPB1_CANFA;ESPB1_CANFA;ESPB1_HUMAN;ESPB1_HUMAN;ESPB1_HUMAN;ESPB1_HUMAN;ESPB1_PIG;ESPB1_PIG;ESPB1_PIG;ESPB1_PIG;FA12_BOVIN;FA12_CAVPO;FA12_HUMAN;FA12_MOUSE;FA12_PIG;FA12_RAT;FINC_BOVIN;FINC_BOVIN;FINC_HUMAN;FINC_HUMAN;FINC_MOUSE;FINC_MOUSE;FINC_NOTVI;FINC_RAT;FINC_RAT;FINC_XENLA;FINC_XENLA;HGFA_CANFA;HGFA_HUMAN;HGFA_MOUSE;LY75_HUMAN;LY75_MESAU;LY75_MOUSE;MMP2_BOVIN;MMP2_BOVIN;MMP2_BOVIN;MMP2_CHICK;MMP2_CHICK;MMP2_CHICK;MMP2_HUMAN;MMP2_HUMAN;MMP2_HUMAN;MMP2_MOUSE;MMP2_MOUSE;MMP2_MOUSE;MMP2_RABIT;MMP2_RABIT;MMP2_RABIT;MMP2_RAT;MMP2_RAT;MMP2_RAT;MMP9_BOVIN;MMP9_BOVIN;MMP9_BOVIN;MMP9_CANFA;MMP9_CANFA;MMP9_CANFA;MMP9_HUMAN;MMP9_HUMAN;MMP9_HUMAN;MMP9_MOUSE;MMP9_MOUSE;MMP9_MOUSE;MMP9_RABIT;MMP9_RABIT;MMP9_RABIT;MMP9_RAT;MMP9_RAT;MMP9_RAT;MPRI_BOVIN;MPRI_HUMAN;MPRI_MOUSE;MRC1L_HUMAN;MRC1_HUMAN;MRC1_MOUSE;MRC2_HUMAN;MRC2_MOUSE;MRC2_RAT;PB1_PIG;PB1_PIG;PLA2R_BOVIN;PLA2R_HUMAN;PLA2R_MOUSE;PLA2R_PONAB;PLA2R_RABIT;SE1L1_HUMAN;SE1L1_MESAU;SE1L1_MOUSE;SE1L1_RAT;SFP1_BOVIN;SFP1_BOVIN;SFP3_BOVIN;SFP3_BOVIN;SFP4_BOVIN;SFP4_BOVIN;SP1_HORSE;SP1_HORSE";
}

#Test parameter based calling
my ($dbobj2,$aln2);
SKIP: {
	ok $dbobj2=Bio::DB::Align::Prosite->new(), 'Bio::DB::Align::Prosite';
	eval{$aln2=$dbobj2->get_Aln_by_acc(-accession=>"PS00023",-format=>"clustalw")};
	skip("Bio::DB::Align::Prosite HTTP error:$@", 7) if $@;
	is $aln2->accession,"PS00023";
	is $aln2->length,"42";
	is $aln2->num_sequences,"87";
	is $aln2->gap_char,"-";
	is $aln2->source,"clustalw";
	
	my $seq=$aln2->get_Seq_by_pos(1);
	is $seq->seq,"CvfPFwYrrliyweCtddgeafgkkWCsltkNFnkdriWkYC";
	is $seq->get_nse,"BSPH1_HUMAN/90-131";
}
