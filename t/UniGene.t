# -*-Perl-*-
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

    $NUMTESTS = 38;
    plan tests => $NUMTESTS;

    unless( eval "require Parse::RecDescent; 1;" ) {
      #warn $@;
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

ok ( defined ($unigene = $str->next_cluster()));

ok($unigene->unigene_id, 'Hs.2');
ok($unigene->title, 'N-acetyltransferase 2 (arylamine N-acetyltransferase)');
ok($unigene->gene, 'NAT2');
ok($unigene->cytoband,'8p22');
ok($unigene->locuslink,'10');
ok($unigene->gnm_terminus,'T');
ok($unigene->chromosome,'8');
ok($unigene->scount,29);
ok(scalar @{ $unigene->express }, 4);
ok(scalar @{ $unigene->sts }, 4);
ok(scalar @{ $unigene->txmap }, 1);
ok(scalar @{ $unigene->protsim } , 4);

ok(scalar @{ $unigene->sequence },29);

ok($unigene->next_express, 'colon');
ok($unigene->next_sts, 'ACC=- NAME=GDB:386004 UNISTS=157141');
ok($unigene->next_txmap, 'D8S549-D8S258; MARKER=stSG40; RHPANEL=GB4');
ok($unigene->next_protsim, 'ORG=Escherischia coli; PROTGI=8928262; PROTID=Ec_pid; PCT=24; ALN=254');
my ($seq1) = $unigene->next_seq;
ok($seq1->display_id, 'D90042');
ok($seq1->desc, 'ACC=D90042 NID=g219415 PID=g219416 UNIGENE_ID=Hs.2');

$unigene->unigene_id('Hs.50');
ok($unigene->unigene_id, 'Hs.50', 'unigene_id was ' . $unigene->unigene_id);

$unigene->title('title_test');
ok($unigene->title, 'title_test', 'title was ' . $unigene->title);

$unigene->gene('gene_test');
ok($unigene->gene, 'gene_test', 'gene was ' . $unigene->gene);

$unigene->cytoband('cytoband_test');
ok($unigene->cytoband, 'cytoband_test', 'cytoband was ' . $unigene->cytoband);

$unigene->locuslink('locuslink_test');
ok($unigene->locuslink, 'locuslink_test', 'locuslink was ' . $unigene->locuslink);

$unigene->gnm_terminus('gnm_terminus_test');
ok($unigene->gnm_terminus, 'gnm_terminus_test', 'gnm_terminus was ' . $unigene->gnm_terminus);

$unigene->chromosome('chromosome_test');
ok($unigene->chromosome, 'chromosome_test', 'chromosome was ' . $unigene->chromosome);

$unigene->scount('scount_test');
ok($unigene->scount, 'scount_test', 'scount was ' . $unigene->scount);

my $seq = $unigene->next_seq;
ok($seq->isa('Bio::PrimarySeqI'), 1,'expected a Bio::PrimarySeq object but got a ' . ref($seq));
my $accession = $seq->accession_number;
ok($accession, 'D90040');

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



