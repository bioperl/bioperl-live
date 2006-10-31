# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;

BEGIN {
	eval { require Test; };
	if ( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 56;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::Annotation::Collection;

ok(1);

# Set to -1 for release version, so warnings aren't printed
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $ast = Bio::SeqIO->new( -format => 'embl',
			   -verbose => $verbose,
			   -file => Bio::Root::IO->catfile
			   ("t","data","roa1.dat"));
$ast->verbose($verbose);
my $as = $ast->next_seq();
ok defined $as->seq;
ok($as->display_id, 'HSHNCPA1');
ok($as->accession_number, 'X79536');
ok($as->seq_version, 1);
ok($as->version, 1);
ok($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
ok($as->molecule, 'RNA');
ok($as->alphabet, 'rna');
ok(scalar $as->all_SeqFeatures(), 4);
ok($as->length, 1198);
ok($as->species->binomial(), 'Homo sapiens');

# EMBL Release 87 changes (8-17-06)

$ast = Bio::SeqIO->new( -format => 'embl',
			   -verbose => $verbose,
			   -file => Bio::Root::IO->catfile
			   ("t","data","roa1_v2.dat"));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok defined $as->seq;
# accession # same as display name now
ok($as->display_id, 'X79536'); 
ok($as->accession_number, 'X79536');
ok($as->seq_version, 1);
ok($as->version, 1);
ok($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
# mRNA instead of RNA
ok($as->molecule, 'mRNA');
ok($as->alphabet, 'rna');
ok(scalar $as->all_SeqFeatures(), 4);
ok($as->length, 1198);
ok($as->species->binomial(), 'Homo sapiens');

my $ent = Bio::SeqIO->new( -file => Bio::Root::IO->catfile
									("t","data","test.embl"),
									-format => 'embl');
my $seq = $ent->next_seq();

ok(defined $seq->seq(), 1,
   'failure to read Embl with ^ location and badly split double quotes');
ok(scalar $seq->annotation->get_Annotations('reference'), 3);

my $out = Bio::SeqIO->new(-file=> ">embl.out",
							  -format => 'embl');
ok($out->write_seq($seq),1,
   'failure to write Embl format with ^ < and > locations');
unlink("embl.out");

# embl with no FT
$ent = Bio::SeqIO->new( -file => Bio::Root::IO->catfile
								("t","data","test.embl"),
								-format => 'embl');


# embl with no FH
my $noFH = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
								("t","data","no_FH.embl"),
								-format => 'embl');
ok(scalar($noFH->next_seq->get_SeqFeatures), 4);


$seq = $ent->next_seq();
ok($seq);
ok(lc($seq->subseq(1,10)),'gatcagtaga');
ok($seq->length);

# bug 1571
$ent = Bio::SeqIO->new(-format => 'embl',
							  -file   => Bio::Root::IO->catfile
							  (qw(t data test.embl2sq)));
ok($ent->next_seq->length,4877);

# embl repbase
$ent = Bio::SeqIO->new(-file => Bio::Root::IO->catfile
							  ("t","data","BEL16-LTR_AG.embl"), -format => 'embl');
$seq = $ent->next_seq;
ok($seq->display_id,'BEL16-LTR_AG');

# test secondary accessions in EMBL (bug #1332)
my $seqio = new Bio::SeqIO(-format => 'embl',
									-file => Bio::Root::IO->catfile
									( qw(t data ECAPAH02.embl)));
$seq = $seqio->next_seq;
ok($seq->accession_number, 'D10483');
ok($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
ok($accs[0], 'J01597');
ok($accs[-1], 'X56742');

### TPA TESTS - Thanks to Richard Adams ###
# test Third Party Annotation entries in EMBL/Gb format 
# to ensure compatability with parsers.
my $str = new Bio::SeqIO(-format =>'embl',
								 -file => Bio::Root::IO->catfile
								 ( qw(t data BN000066-tpa.embl)));
$seq = $str->next_seq;
ok(defined $seq);
ok($seq->accession_number, 'BN000066');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'AGA000066');
ok($seq->length, 5195);
ok($seq->division, 'INV');
ok($seq->get_dates, 2);
ok($seq->keywords, 'acetylcholinesterase; achE1 gene; Third Party Annotation; TPA.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 15);

my $spec_obj = $seq->species;
ok ($spec_obj->common_name, 'African malaria mosquito');
ok ($spec_obj->species, 'gambiae');
ok ($spec_obj->genus, 'Anopheles');
ok ($spec_obj->binomial, 'Anopheles gambiae');

my $ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[1];
ok ($reference->title,'"A novel acetylcholinesterase gene in mosquitoes codes for the insecticide target and is non-homologous to the ace gene in Drosophila"');
ok ($reference->authors,'Weill M., Fort P., Berthomi eu A., Dubois M.P., Pasteur N., Raymond M.');
my $cmmnt =  ($ac->get_Annotations('comment') )[0];
ok($cmmnt->text, 'see also AJ488492 for achE-1 from Kisumu strain Third Party Annotation Database: This TPA record uses Anopheles gambiae trace archive data (http://trace.ensembl.org) ');


$ent = Bio::SeqIO->new( -file => Bio::Root::IO->catfile
                        ("t","data","test.embl"),
                        -format => 'embl');
$ent->verbose($verbose);
$seq = $ent->next_seq();
my $species = $seq->species();
my @cl = $species->classification();
ok( $cl[3] ne $species->genus(), 1, 'genus duplicated in EMBL parsing');
$ent->close();

#
## read-write - test embl writing of a PrimarySeq
#
my $primaryseq = new Bio::PrimarySeq( -seq => 'AGAGAGAGATA',
                                      -id  => 'myid',
                                      -desc => 'mydescr',
                                      -alphabet => 'DNA',
                                      -accession_number => 'myaccession');



$verbose = -1 unless $ENV{'BIOPERLDEBUG'};  # silence warnings unless we are debuggin

my $embl = new Bio::SeqIO(-format => 'embl',
                          -verbose => $verbose,
                          -file => ">primaryseq.embl");

ok($embl->write_seq($primaryseq));

# this should generate a warning
my $scalar = "test";
eval {
	$embl->write_seq($scalar);
};
ok ($@);

unlink("primaryseq.embl");
