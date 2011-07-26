# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

#test Bio::DB::Align::Pfam
BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 21,
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
	ok $dbobj=Bio::DB::Align->new(-db=>"Rfam"), 'Bio::DB::Align::Rfam';
	eval{$aln=$dbobj->get_Aln_by_acc("RF00360")};
	skip("Bio::DB::Align::Rfam HTTP error:$@", 5) if $@;
	is $aln->accession,"RF00360";
	is $aln->length,"132";
	is $aln->num_sequences,"9";
	is $aln->source,"Rfam";
	
	my @names;
	foreach my $seq ($aln->next_Seq) {
		push @names,$seq->display_id;
	}
	is join(";",@names),"O.sativa.1;O.sativa.2;O.sativa.3;O.sativa.4;O.sativa.5;O.sativa.6;M.truncatula.1;A.thaliana.1;A.thaliana.2";
}

#Test parameter based calling
my ($dbobj2,$aln2);
SKIP: {
	ok $dbobj2=Bio::DB::Align::Rfam->new(), 'Bio::DB::Align::Rfam';
	eval{$aln2=$dbobj2->get_Aln_by_acc(-accession=>"RF00360",-alignment=>"full",-nselabel=>1, -format=>"stockholm",-gap=>"dashes")};
	skip("Bio::DB::Align::Rfam HTTP error:$@", 7) if $@;
	is $aln2->accession,"RF00360";
	is $aln2->length,"162";
	is $aln2->num_sequences,"78";
	is $aln2->gap_char,"-";
	is $aln2->source,"stockholm";
	
	my $seq=$aln2->get_Seq_by_pos(1);
	is $seq->seq,"GU-UUGCAGUGACGACAAGAAAAUUUCGUCAAGCUCAACAGACUUGAUUACGGGGAUAG--------------GAACAUCGUUUCGCAUCCCAUAUGAUAGUUAACCCGCUGAUCCGA--------------------------------------GCAAAC";
	is $seq->display_id,"ABEU01007475";
}

#Test nongapped retrieval
my $aln3;
SKIP: {
	eval{$aln3=$dbobj->get_Aln_by_acc(-accession=>"RF01050", -gap=>"none")};
	skip("Bio::DB::Align::Rfam HTTP error:$@", 5) if $@;
	is $aln3->accession,"RF01050";
	is $aln3->length,"1220";
	is $aln3->num_sequences,"13"; 
	my $seq2=$aln3->get_Seq_by_pos(1);
	is $seq2->seq,"CAAAGGAAGAUAGGUAACCUAUUAAGAUGUCAGCGGCUGUUGCGUUUGCUUAGUUGUUUUUUUUUUUAGUAUUUGUUUUUUGUACAUUUUCCGUUUGAAUUUUCUAUCAUGCAAGCCUCAGAGAUUUGGUAGAUGCUUAAUGAUGUAAAGGUUGCGUUGAAUUUUGAAUUGUUCUCUCAAGUUAGCAGGCAGGUGCACUUUUCUCUUUUACGAAGGACUUCGAGUGCACUGGGGUCAGCAUAGAUACUUGUGUGUAUAUAUAUUUGUGGUUUUUGUUAUUUUUCUACUUAUAGAUGGCUAAAAUCUGAGUUUGAAAGAUGGCCACCAUAAAUUCUUAAAAGUGGUAUCGCAUUUAGCCCCCCUUCGUACCAAUUCUGUUUUCGUCUCAAGCUCUUCAUUACUACAGCGCAAGUCUACCAUUACCACACCCACACACAGAUAUCACGGCUAAUGUUUAUUAGUUAAGUUUCCAUGAGCUCACUCUUUAUUCUUUUCCUCGUUUUCUUAUACCUAGUAUGUUUUCUGACGCUUUUUGAAGUGACAGAAAAAAGGAGUUUAAAUUAGAUUUGCAAACGGACGGUACUAAAUACCGUCACUUUUACGUCUAACUUAUCGUUAACUCUGCAAAAAGAAAACAAGAAAAAGAAAAUCAGUGAAUAGGAGUAUAUAUAGAAAUGGUUUAUUCUAUUUUUUUUUCGUUUUUCAGUAGAUUUUCGCCUUUAAAAGAAUAAAUUCCAUUAAAAAAAGGUAAAAAGAAAAAUCUAUUCACUGAACUUGCCGAUAGAAUUUGCAAAUGUGUCAAGUGCAUCAAGGAGAAGAUAGAGGAGAAACUCGAUUAGGCAAACAAGCCAAAAGGCAAGGAUACCCCUUCUCAAGCAUUGUUAGGUUUAUGGGCUACCAGUAAGCGAAAAAUGAUACAGGACGAAAAAAACAUAAUUCGAGAUUUUUCAGGAUGGAUUUUUUUAGGUAUCUAUCAAAACUAUUCUGAUGAUCAAUACAGUAUUUUUGUCGCAUUAUUAUUACGGUAGUGUAGACAUUAGUUUCCAAGCGGAAGAAACUGUGUGUUCAUUUUAUGGAUUUUCGUGUUGUACGAUUUUUUUCAGUUGCGUUAGCAACUACGUGCCCACAAUACUUCCGAUGCAUUUAGAUAAUUUUUGGAAACAUU-------------------------------------------------------";
	is $seq2->get_nse,"S.paradoxus.1/1-1165";
}
