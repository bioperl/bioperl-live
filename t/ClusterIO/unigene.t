# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 73);
	
	use_ok('Bio::ClusterIO');
}

my ($str, $unigene); # predeclare variables for strict

$str = Bio::ClusterIO->new('-file' => test_input_file('unigene.data'), '-format' => "unigene");
ok $str, 'new Bio::ClusterIO object defined';

ok ( defined ($unigene = $str->next_cluster()));

# check interface implementations to be sure
isa_ok $unigene, "Bio::Cluster::UniGeneI";
isa_ok $unigene, "Bio::ClusterI";
isa_ok $unigene, "Bio::IdentifiableI";
isa_ok $unigene, "Bio::DescribableI";

# test specific instance of unigene record provided in the unigene.data file
is ($unigene->unigene_id, 'Hs.2');
is ($unigene->title, 'N-acetyltransferase 2 (arylamine N-acetyltransferase)');
is ($unigene->gene, 'NAT2');
is ($unigene->cytoband,'8p22');
is ($unigene->gnm_terminus,'S');
is ($unigene->homol,'YES');
is ($unigene->restr_expr,'liver');
is ($unigene->scount,26);
is (scalar @{ $unigene->locuslink }, 1);
is (scalar @{ $unigene->chromosome }, 1);
is (scalar @{ $unigene->express }, 7);
is (scalar @{ $unigene->sts }, 8);
is (scalar @{ $unigene->txmap }, 0);
is (scalar @{ $unigene->protsim } , 4);
is (scalar @{ $unigene->sequences },26);

is ($unigene->locuslink->[0], '10');
is ($unigene->chromosome->[0], '8');
is ($unigene->express->[0], 'liver');
is ($unigene->sts->[0], 'ACC=G59899 UNISTS=137181');
is ($unigene->protsim->[0], 'ORG=Escherischia coli; PROTGI=16129422; PROTID=ref:NP_415980.1; PCT=24.81; ALN=255');

my ($seq1) = $unigene->next_seq;
is ($seq1->display_id, 'BX095770');

# test recognition of species

SKIP: {
	# TODO? why wouldn't it be available? Bug if not? Remove this skip?
	skip 'species not available', 'species test' unless defined $unigene->species;
	is $unigene->species->binomial, "Homo sapiens" ;
}


# test accessors of interfaces
is  ($seq1->namespace, "GenBank");
is  ($seq1->authority, "NCBI");
is  ($seq1->alphabet, "dna");
my $n = 1; # we've seen already one seq
while($seq1 = $unigene->next_seq()) {
    $n++;
}
is  ($n, 26);
is  ($unigene->size(), 26);
is  (scalar($unigene->get_members()), 26);
is  ($unigene->description, 'N-acetyltransferase 2 (arylamine N-acetyltransferase)');
is  ($unigene->display_id, "Hs.2");
is  ($unigene->namespace, "UniGene");
is  ($unigene->authority, "NCBI");

$unigene->unigene_id('Hs.50');
is ($unigene->unigene_id, 'Hs.50') || diag('unigene_id was ' . $unigene->unigene_id);

$unigene->title('title_test');
is ($unigene->title, 'title_test') || diag('title was ' . $unigene->title);

$unigene->gene('gene_test');
is ($unigene->gene, 'gene_test') || diag('gene was ' . $unigene->gene);

$unigene->cytoband('cytoband_test');
is ($unigene->cytoband, 'cytoband_test') || diag('cytoband was ' . $unigene->cytoband);

$unigene->gnm_terminus('gnm_terminus_test');
is ($unigene->gnm_terminus, 'gnm_terminus_test') || diag('gnm_terminus was ' . $unigene->gnm_terminus);

$unigene->homol('homol_test');
is ($unigene->homol, 'homol_test') || diag('homol was ' . $unigene->homol);

$unigene->restr_expr('restr_expr_test');
is ($unigene->restr_expr, 'restr_expr_test') || diag('restr_expr was ' . $unigene->restr_expr);

$unigene->scount('scount_test');
is ($unigene->scount, 'scount_test') || diag('scount was ' . $unigene->scount);

my $seq = $unigene->next_seq;
$seq = $unigene->next_seq;
isa_ok ($seq, 'Bio::PrimarySeqI') || diag('expected a Bio::PrimarySeq object but got a ' . ref($seq));
my $accession = $seq->accession_number;
is ($accession, 'AI262683');
my $version = $seq->seq_version();
is ($version, 1);

# test the sequence parsing is working
my $ac = $seq->annotation();
my $simple_ann_object;
($simple_ann_object) = $ac->get_Annotations('seqtype');
ok defined $simple_ann_object, 'annotation object defined';
is ($simple_ann_object->value(), 'EST') || diag('seqtype was ' . $simple_ann_object->value);	

# test PERIPHERAL, bug 1708
$seq = $unigene->next_seq;
$accession = $seq->accession_number;
is ($accession, 'CB161982');

my @acs = $seq->annotation->get_Annotations('peripheral');
is  $acs[0]->display_text, 1;

# tests not specific to unigene record provided in the unigene.data file
my @locuslink_test = ( "58473", "5354" );
$unigene->locuslink(\@locuslink_test);
my @locuslink_results;
while (my $locuslink = $unigene->next_locuslink) {
	push @locuslink_results, $locuslink;
}
is(scalar(@locuslink_results), 2) || diag('expected locuslink to have 2 entries but it had ' . scalar(@locuslink_results));
my $locuslink = shift @locuslink_results;
is( $locuslink, '58473') || diag('expected 58473 but got ' . $locuslink);


my @express_test = qw( kidney heart liver spleen );
$unigene->express(\@express_test);
my @express_results;
while (my $tissue = $unigene->next_express) {
	push @express_results, $tissue;
}
is(  scalar(@express_results), 4) || diag('expected express to have 4 entries but it had ' . scalar(@express_results));

my @chromosome_test = ( "7", "11" );
$unigene->chromosome(\@chromosome_test);
my @chromosome_results;
while (my $chromosome = $unigene->next_chromosome) {
	push @chromosome_results, $chromosome;
}
is( scalar(@chromosome_results), 2) || diag('expected chromosome to have 2 entries but it had ' . scalar(@chromosome_results));
my $chromosome = shift @chromosome_results;
is(  $chromosome, '7') || diag('expected 7 but got ' . $chromosome);

my @sts_test = ( "ACC=- NAME=sts-D90276 UNISTS=37687", "ACC=G29786 NAME=SHGC-35230 UNISTS=58455" );
$unigene->sts(\@sts_test);
my @sts_results;
while (my $sts = $unigene->next_sts) {
	push @sts_results, $sts;
}
is(scalar(@sts_results), 2) || diag('expected sts to have 2 entries but it had ' . scalar(@sts_results));
my $sts = shift @sts_results;
is($sts, 'ACC=- NAME=sts-D90276 UNISTS=37687') || diag('expected ACC=- NAME=sts-D90276 UNISTS=37687 but got ' . $sts);

my @txmap_test = ("D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4" , "D19S425-D19S418; MARKER=stSG41396; RHPANEL=GB4");
$unigene->txmap(\@txmap_test);
my @txmap_results;
while (my $txmap = $unigene->next_txmap) {
	push @txmap_results, $txmap;
}
is(scalar(@txmap_results), 2) || diag('expected txmap to have 2 entries but it had ' . scalar(@txmap_results));
my $txmap = shift @txmap_results;
is ($txmap, 'D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4') || diag('expected D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4 but got ' . $txmap);

my @protsim_test = ("ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243" , "ORG=Mus musculus; PROTGI=2497288; PROTID=sp:Q61400; PCT=42; ALN=143");
$unigene->protsim(\@protsim_test);
my @protsim_results;
while (my $protsim = $unigene->next_protsim) {
    push @protsim_results, $protsim;
}
is (scalar(@protsim_results), 2) || diag('expected protsim to have 2 entries but it had ' . scalar(@protsim_results));
my $protsim = shift @protsim_results;
is ($protsim, 'ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243') || diag('expected ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243 but got ' . $protsim);



# do a quick test on Rn record included as the next cluster in the
# test data file because it has version numbers tacked on the end of
# the accession numbers in each seq line - NCBI has started doing this
# now (Sept 2003).

$unigene = $str->next_cluster();
$seq = $unigene->next_seq;
isa_ok ($seq,'Bio::PrimarySeqI') || diag( 'expected a Bio::PrimarySeq object but got a ' . ref($seq));
$version = $seq->seq_version();
is($version, '1');

# next cluster contains a // in the title - yes NCBI did that. Nonetheless,
# this should not trip up the parser:

$unigene = $str->next_cluster();
ok ($unigene, 'next cluster'); # previously this would have been undef
is  ($unigene->unigene_id, "Mm.340763");
is ($unigene->title, 'Transcribed locus, strongly similar to NP_003008.1 splicing factor, arginine/serine-rich 3; splicing factor, arginine//serine-rich, 20-kD [Homo sapiens]');
is ($unigene->homol, 'YES');
is ($unigene->scount, 31);
is (scalar($unigene->get_members()), 31);
