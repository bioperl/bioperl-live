## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
use lib '.','./blib/lib';

my $error;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 37;
    plan tests => $NUMTESTS;

    unless( eval "use Parse::RecDescent; 1;" ) {
      warn $@;
      print STDERR "Parse::RecDescent not installed. This means that Bio::Cluster::UniGene module is not usable. Skipping tests.\n";
      for( 1..$NUMTESTS ) {
	skip(1,"Parse::RecDescent not installed. This means that Bio::Cluster::UniGene module is not usable. Skipping tests.\n");
      }
      $error = 1;
    }
}

if( $error ==  1 ) {
    exit(0);
}

use Bio::Cluster::UniGene;
use Bio::ClusterIO;

my ($str, $unigene); # predeclare variables for strict



$str = Bio::ClusterIO->new('-file' => Bio::Root::IO->catfile("t","data","unigene.data"), '-format' => "unigene");

ok $str;

ok ( defined ($unigene = $str->next_unigene()));

ok $unigene->unigene_id;
ok $unigene->title;
ok $unigene->gene;
ok $unigene->cytoband;
ok $unigene->locuslink;
ok $unigene->gnm_terminus;
ok $unigene->chromosome;
ok $unigene->scount;
ok $unigene->express;
ok $unigene->sts;
ok $unigene->txmap;
ok $unigene->protsim;
ok $unigene->sequence;

ok $unigene->next_express;
ok $unigene->next_sts;
ok $unigene->next_txmap;
ok $unigene->next_protsim;
ok $unigene->next_seq;

$unigene->unigene_id('Hs.50');
ok $unigene->unigene_id, 'Hs.50', 'unigene_id was ' . $unigene->unigene_id;

$unigene->title('title_test');
ok $unigene->title, 'title_test', 'title was ' . $unigene->title;

$unigene->gene('gene_test');
ok $unigene->gene, 'gene_test', 'gene was ' . $unigene->gene;

$unigene->cytoband('cytoband_test');
ok $unigene->cytoband, 'cytoband_test', 'cytoband was ' . $unigene->cytoband;

$unigene->locuslink('locuslink_test');
ok $unigene->locuslink, 'locuslink_test', 'locuslink was ' . $unigene->locuslink;

$unigene->gnm_terminus('gnm_terminus_test');
ok $unigene->gnm_terminus, 'gnm_terminus_test', 'gnm_terminus was ' . $unigene->gnm_terminus;

$unigene->chromosome('chromosome_test');
ok $unigene->chromosome, 'chromosome_test', 'chromosome was ' . $unigene->chromosome;

$unigene->scount('scount_test');
ok $unigene->scount, 'scount_test', 'scount was ' . $unigene->scount;



my $seq = $unigene->next_seq;
ok ref($seq),qr/Bio::PrimarySeq/, 'expected a Bio::PrimarySeq object but got a ' . ref($seq);
my $accession = $seq->accession_number;
ok $accession;

my @express_test = qw( kidney heart liver spleen );
$unigene->express(\@express_test);
my @express_results;
while (my $tissue = $unigene->next_express) {
	push @express_results, $tissue;
}
ok scalar(@express_results), 4, 'expected express to have 4 entries but it had ' . scalar(@express_results);



my @sts_test = ( "ACC=- NAME=sts-D90276 UNISTS=37687", "ACC=G29786 NAME=SHGC-35230 UNISTS=58455" );
$unigene->sts(\@sts_test);
my @sts_results;
while (my $sts = $unigene->next_sts) {
	push @sts_results, $sts;
}
ok scalar(@sts_results), 2, 'expected sts to have 2 entries but it had ' . scalar(@sts_results);
my $sts = shift @sts_results;
ok $sts, 'ACC=- NAME=sts-D90276 UNISTS=37687', 'expected ACC=- NAME=sts-D90276 UNISTS=37687 but got ' . $sts;



my @txmap_test = ("D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4" , "D19S425-D19S418; MARKER=stSG41396; RHPANEL=GB4");
$unigene->txmap(\@txmap_test);
my @txmap_results;
while (my $txmap = $unigene->next_txmap) {
	push @txmap_results, $txmap;
}
ok scalar(@txmap_results), 2, 'expected txmap to have 2 entries but it had ' . scalar(@txmap_results);
my $txmap = shift @txmap_results;
ok $txmap, 'D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4', 'expected D19S425-D19S418; MARKER=sts-D90276; RHPANEL=GB4 but got ' . $txmap;



my @protsim_test = ("ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243" , "ORG=Mus musculus; PROTGI=2497288; PROTID=sp:Q61400; PCT=42; ALN=143");
$unigene->protsim(\@protsim_test);
my @protsim_results;
while (my $protsim = $unigene->next_protsim) {
	push @protsim_results, $protsim;
}
ok scalar(@protsim_results), 2, 'expected protsim to have 2 entries but it had ' . scalar(@protsim_results);
my $protsim = shift @protsim_results;
ok $protsim, 'ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243', 'expected ORG=Homo sapiens; PROTGI=107211; PROTID=pir:A40428; PCT=100; ALN=243 but got ' . $protsim;



