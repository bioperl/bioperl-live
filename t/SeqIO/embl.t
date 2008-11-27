# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
	use BioperlTest;
	
	test_begin(-tests => 56);
	
    use_ok('Bio::SeqIO::embl');
}

my $verbose = test_debug();

my $ast = Bio::SeqIO->new( -format => 'embl',
			   -verbose => $verbose,
			   -file => test_input_file('roa1.dat'));
$ast->verbose($verbose);
my $as = $ast->next_seq();
ok defined $as->seq;
is($as->display_id, 'HSHNCPA1');
is($as->accession_number, 'X79536');
is($as->seq_version, 1);
is($as->version, 1);
is($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
is($as->molecule, 'RNA');
is($as->alphabet, 'rna');
is(scalar $as->all_SeqFeatures(), 4);
is($as->length, 1198);
is($as->species->binomial(), 'Homo sapiens');

# EMBL Release 87 changes (8-17-06)

$ast = Bio::SeqIO->new( -format => 'embl',
			   -verbose => $verbose,
			   -file => test_input_file('roa1_v2.dat'));
$ast->verbose($verbose);
$as = $ast->next_seq();
ok defined $as->seq;
# accession # same as display name now
is($as->display_id, 'X79536'); 
is($as->accession_number, 'X79536');
is($as->seq_version, 1);
is($as->version, 1);
is($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
# mRNA instead of RNA
is($as->molecule, 'mRNA');
is($as->alphabet, 'rna');
is(scalar $as->all_SeqFeatures(), 4);
is($as->length, 1198);
is($as->species->binomial(), 'Homo sapiens');

my $ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
			   -format => 'embl');
my $seq = $ent->next_seq();

is(defined $seq->seq(), 1,
   'success reading Embl with ^ location and badly split double quotes');
is(scalar $seq->annotation->get_Annotations('reference'), 3);

my $out_file = test_output_file();
my $out = Bio::SeqIO->new(-file=> ">$out_file",
			  -format => 'embl');
is($out->write_seq($seq),1,
   'success writing Embl format with ^ < and > locations');

# embl with no FT
$ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
			-format => 'embl');
$seq = $ent->next_seq();

ok($seq);
is(lc($seq->subseq(1,10)),'gatcagtaga');
is($seq->length, 4870);

# embl with no FH
my $noFH = Bio::SeqIO->new(-file => test_input_file('no_FH.embl'),
			-format => 'embl');
is(scalar($noFH->next_seq->get_SeqFeatures), 4);

# bug 1571
$ent = Bio::SeqIO->new(-format => 'embl',
		       -file   => test_input_file('test.embl2sq'));
is($ent->next_seq->length,4877);

# embl repbase
$ent = Bio::SeqIO->new(-file => test_input_file('BEL16-LTR_AG.embl'), -format => 'embl');
$seq = $ent->next_seq;
is($seq->display_id,'BEL16-LTR_AG');

# test secondary accessions in EMBL (bug #1332)
my $seqio = Bio::SeqIO->new(-format => 'embl',
			   -file => test_input_file('ECAPAH02.embl'));
$seq = $seqio->next_seq;
is($seq->accession_number, 'D10483');
is($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
is($accs[0], 'J01597');
is($accs[-1], 'X56742');

### TPA TESTS - Thanks to Richard Adams ###
# test Third Party Annotation entries in EMBL/Gb format 
# to ensure compatability with parsers.
my $str = Bio::SeqIO->new(-format =>'embl',
			 -file => test_input_file('BN000066-tpa.embl'));
$seq = $str->next_seq;
ok(defined $seq);
is($seq->accession_number, 'BN000066');
is($seq->alphabet, 'dna');
is($seq->display_id, 'AGA000066');
is($seq->length, 5195);
is($seq->division, 'INV');
is($seq->get_dates, 2);
is($seq->keywords, 'acetylcholinesterase; achE1 gene; Third Party Annotation; TPA.');
is($seq->seq_version, 1);
is($seq->feature_count, 15);

my $spec_obj = $seq->species;
is ($spec_obj->common_name, 'African malaria mosquito');
is ($spec_obj->species, 'gambiae');
is ($spec_obj->genus, 'Anopheles');
is ($spec_obj->binomial, 'Anopheles gambiae');

my $ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[1];
is ($reference->title,'"A novel acetylcholinesterase gene in mosquitoes codes for the insecticide target and is non-homologous to the ace gene in Drosophila"');
is ($reference->authors,'Weill M., Fort P., Berthomi eu A., Dubois M.P., Pasteur N., Raymond M.');
my $cmmnt =  ($ac->get_Annotations('comment') )[0];
is($cmmnt->text, 'see also AJ488492 for achE-1 from Kisumu strain Third Party Annotation Database: This TPA record uses Anopheles gambiae trace archive data (http://trace.ensembl.org) ');


$ent = Bio::SeqIO->new( -file => test_input_file('test.embl'),
                        -format => 'embl');
$ent->verbose($verbose);
$seq = $ent->next_seq();
my $species = $seq->species();
my @cl = $species->classification();
is( $cl[3] ne $species->genus(), 1, 'genus duplication test');
$ent->close();

#
## read-write - test embl writing of a PrimarySeq
#
my $primaryseq = Bio::PrimarySeq->new( -seq => 'AGAGAGAGATA',
                                      -id  => 'myid',
                                      -desc => 'mydescr',
                                      -alphabet => 'DNA',
                                      -accession_number => 'myaccession');



$verbose = -1 unless $ENV{'BIOPERLDEBUG'};  # silence warnings unless we are debuggin

my $embl = Bio::SeqIO->new(-format => 'embl',
                          -verbose => $verbose,
                          -file => ">$out_file");

ok($embl->write_seq($primaryseq));

# this should generate a warning
my $scalar = "test";
eval {
	$embl->write_seq($scalar);
};
ok ($@);
