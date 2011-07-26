# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

#test Bio::DB::Align::Pfam
BEGIN {
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 20,
			   -requires_modules => [qw(HTTP::Request
									    LWP::UserAgent
										Bio::AlignIO
										Bio::DB::GenericWebAgent
										Bio::DB::Align
										Bio::DB::EUtilities)],
			   -requires_networking => 1,
			   -requires_email      => 1,);
	
	use_ok('Bio::Root::IO');
	use_ok('Bio::Root::Root');
}

my $email = test_email();

#test calling from Bio::DB::Align
SKIP: {
	my ($dbobj, $aln);
	ok $dbobj=Bio::DB::Align->new(-db=>"ProtClustDB",-email=>$email), 'Bio::DB::Align::ProtClustDB';
	eval{$aln=$dbobj->get_Aln_by_id("2725839")};
	skip("Bio::DB::Align::ProtClustDB HTTP error:$@", 8) if $@;
	
	is $aln->id,"2725839";
	is $aln->accession,"CLSN2725839";
	is $aln->length,"124";
	is $aln->num_sequences,"2";
	is $aln->source,"ProtClustDB";
	
	my @names;
	foreach my $seq ($aln->next_Seq) {
		push @names,$seq->display_id;
	}
	is join(";",@names),"ref|XP_002515767;ref|XP_002519568";
	
	my $acc=$dbobj->id2acc("2725839");
	is $acc,"CLSN2725839";
	my $id=$dbobj->acc2id("CLSN2725839");
	is $id,"2725839";
}

#Test parameter based calling
SKIP: {
	my ($dbobj2,$aln2);
	ok $dbobj2=Bio::DB::Align::ProtClustDB->new(-email=>$email);
	eval{$aln2=$dbobj2->get_Aln_by_acc(-accession=>"CLSN2725839")};
	skip("Bio::DB::Align::ProtClustDB HTTP error:$@", 8) if $@;
	is $aln2->id,"2725839";
	is $aln2->accession,"CLSN2725839";
	is $aln2->length,"124";
	is $aln2->num_sequences,"2";
	is $aln2->gap_char,"-";
	is $aln2->source,"ProtClustDB";
	
	my $seq=$aln2->get_Seq_by_pos(1);
	is $seq->seq,"-------------------------------------------------------------------MRRLRGSNVELFKRPQCLESKGHVVFIRADVNQPGSYNSTSPSLAVVIVTMNWPIAN";
	is $seq->get_nse,"ref|XP_002515767/1-57";
}