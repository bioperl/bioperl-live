# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan test => 134;
}

use Bio::SearchIO;
use Bio::Tools::HMMER::Domain;
use Bio::Tools::HMMER::Set;
use Bio::Tools::HMMER::Results;
use Bio::Root::IO;


my $searchio = new Bio::SearchIO(-format => 'hmmer',
				 -file   => Bio::Root::IO->catfile
				 ("t","data","hmmpfam.out"));

while( my $result = $searchio->next_result ) {
    ok(ref($result),'Bio::Search::Result::HMMERResult');
    ok($result->algorithm, 'HMMPFAM');
    ok($result->algorithm_version, '2.1.1');
    ok($result->hmm_name, 'pfam');
    ok($result->sequence_file, '/home/birney/src/wise2/example/road.pep');
    ok($result->query_name, 'roa1_drome');
    ok($result->query_description, '');
    ok($result->num_hits(), 1);
    while( my $hit = $result->next_model ) {
	ok($hit->name, 'SEED');
	ok($hit->raw_score, '146.1');
	ok($hit->significance, '6.3e-40');
	ok(ref($hit), 'Bio::Search::Hit::HMMERHit');
	ok($hit->num_hsps, 2);
	my $hsp = $hit->next_domain;
	if( defined $hsp ) {
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 77);
	    ok($hsp->query->start, 33);
	    ok($hsp->query->end, 103);
	    ok($hsp->score, 71.2);
	    ok($hsp->evalue, '2.2e-17');
	    ok($hsp->query_string, 'LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP');
	    ok($hsp->gaps('query'), 7);
	    ok($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv');
	    ok($hsp->homology_string, 'lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ');
	    ok(	length($hsp->homology_string), length($hsp->hit_string));
	    ok( length($hsp->query_string), length($hsp->homology_string));
	}
	$hsp = $hit->next_domain;
	if( defined $hsp ) {
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 77);
	    ok($hsp->query->start, 124);
	    ok($hsp->query->end, 194);
	    ok($hsp->score, 75.5);
	    ok($hsp->evalue, '1.1e-18');
	    ok($hsp->query_string, 'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV');
	    ok($hsp->gaps('query'), 6);
	    ok($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv');	 
	    ok($hsp->homology_string, 'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v');
	    ok(	length($hsp->homology_string), length($hsp->hit_string));
	    ok( length($hsp->query_string), length($hsp->homology_string));
	}
    }
}
$searchio = new Bio::SearchIO(-format => 'hmmer',
			      -file   => Bio::Root::IO->catfile
			      ("t","data","hmmsearch.out"));
while( my $result = $searchio->next_result ) {
    ok(ref($result),'Bio::Search::Result::HMMERResult');
    ok($result->algorithm, 'HMMSEARCH');
    ok($result->algorithm_version, '2.0');
    ok($result->hmm_name, 'HMM [SEED]');
    ok($result->sequence_file, 'HMM.dbtemp.29591');
    ok($result->query_name, 'SEED');
    ok($result->query_description, '');
    ok($result->num_hits(), 1215);
    my $hit = $result->next_model;
    ok($hit->name, 'Q91581');
    ok($hit->description,'Q91581 POLYADENYLATION FACTOR 64 KDA SUBUN');
    ok($hit->significance, '2e-31');
    ok($hit->raw_score, 119.7);
    my $hsp = $hit->next_domain;
    ok($hsp->score,119.7);
    ok($hsp->evalue, '2e-31');
    ok($hsp->query->start, 18);
    ok($hsp->query->end, 89);
    ok($hsp->hit->start, 1);
    ok($hsp->hit->end, 77);
    ok($hsp->query->seqname(), 'SEED');
    ok($hsp->hit->seqname(), 'Q91581');   
}

$searchio = new Bio::SearchIO(-format => 'hmmer',
			      -file   => Bio::Root::IO->catfile("t","data",
								"L77119.hmmer"));

while( my $result = $searchio->next_result ) {
    ok(ref($result),'Bio::Search::Result::HMMERResult');
    ok($result->algorithm, 'HMMPFAM');
    ok($result->algorithm_version, '2.2g');
    ok($result->hmm_name, 'Pfam');
    ok($result->sequence_file, 'L77119.faa');
    ok($result->query_name, 'gi|1522636|gb|AAC37060.1|');
    ok($result->query_description, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]');
    ok($result->num_hits(), 1);
    my $hit = $result->next_hit;
    ok($hit->name, 'Methylase_M');
    ok($hit->description,'Type I restriction modification system, M');
    ok($hit->significance, '0.0022');
    ok($hit->raw_score, -105.2);
    my $hsp = $hit->next_hsp;
    ok($hsp->score,-105.2);
    ok($hsp->evalue, '0.0022');
    ok($hsp->query->start, 280);
    ok($hsp->query->end, 481);
    ok($hsp->hit->start, 1);
    ok($hsp->hit->end, 279);
    ok($hsp->query->seqname(), 'gi|1522636|gb|AAC37060.1|');
    ok($hsp->hit->seqname(), 'Methylase_M');
    ok($hsp->hit_string, 'lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn');
    ok($hsp->query_string, 'NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST');
    ok($hsp->homology_string, ' ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   ++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf ++++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+');
    
}


$searchio = new Bio::SearchIO(-format => 'hmmer',
			      -file   => Bio::Root::IO->catfile("t","data",
								"cysprot1b.hmmsearch"));


while( my $result = $searchio->next_result ) {
    ok(ref($result),'Bio::Search::Result::HMMERResult');
    ok($result->algorithm, 'HMMSEARCH');
    ok($result->algorithm_version, '2.2g');
    ok($result->hmm_name, 'Peptidase_C1.hmm [Peptidase_C1]');
    ok($result->sequence_file, 'cysprot1b.fa');
    ok($result->query_name, 'Peptidase_C1');
    ok($result->query_accession, 'PF00112');
    ok($result->query_description, 'Papain family cysteine protease');
    ok($result->num_hits(), 4);
    my $hit = $result->next_hit;
    ok($hit->name, 'CATL_RAT');
    ok($hit->description,'');
    ok($hit->significance, '2e-135');
    ok($hit->raw_score, 449.4);
    my $hsp = $hit->next_hsp;
    ok($hsp->score,449.4);
    ok($hsp->evalue, '2e-135');
    ok($hsp->query->start, 1);
    ok($hsp->query->end, 337);
    ok($hsp->hit->start, 114);
    ok($hsp->hit->end, 332);
    ok($hsp->query->seqname(), 'Peptidase_C1');
    ok($hsp->hit->seqname(), 'CATL_RAT');
    ok($hsp->hit_string, 'IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQMFLKT------GKLISLSEQNLVDCSH-DQGNQ------GCNG-GLMDFAFQYIKE-----NGGLDSEESY-----PYE----AKD-------------------GSCKYR-AEYAV-----ANDTGFVDIPQQ-----EKALMKAVATVGPISVAMDASHPS---LQFYSSG-------IYYEP---NCSSK---DLDHGVLVVGYGYEG-T------------------------------------DSNKDKYWLVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI');
    ok($hsp->homology_string, '+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+ ++kt      gkl+sLSEQ+LvDC++ d gn+      GCnG Glmd Af+Yik+     NgGl++E++Y     PY+    +kd                   g+Cky+  + ++     a+++g++d+p++     E+al+ka+a++GP+sVa+das+ s    q+Y+sG       +Y+++    C+++   +LdH+Vl+VGYG e+                                      ++++ +YW+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi');
    ok($hsp->query_string, 'lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgtkawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikkeqIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgtCkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVaidasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGYGteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYWIVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi');
    $hit = $result->next_hit;
    ok($hit->name, 'CATL_HUMAN');
    ok($hit->description,'');
    ok($hit->significance, '6.1e-134');
    ok($hit->raw_score, 444.5);
}
			      
my ($domain,$set,$homol,$rev,$res,$dom,@doms);
    $domain = Bio::Tools::HMMER::Domain->new(-verbose=>1);

ok ref($domain), 'Bio::Tools::HMMER::Domain';

$domain->start(50);
$domain->end(200);
$domain->hstart(10);
$domain->hend(100);
$domain->seqbits(50);
$domain->bits(20);
$domain->evalue(0.0001);
$domain->seqname('silly');


# test that we can get out forward and reverse homol_SeqFeatures
$homol = $domain->feature2();
ok $homol->start(), 10;

$rev = $domain;

ok $rev->start(), 50;

$set = Bio::Tools::HMMER::Set->new();
$set->add_Domain($domain);

@doms = $set->each_Domain();
$dom = shift @doms;

ok $dom->start(), 50;

$set->bits(300);
$set->evalue(0.0001);
$set->name('sillyname');
$set->desc('a desc');
$set->accession('fakeaccesssion');
ok $set->bits(), 300;
ok $set->evalue(), 0.0001;
ok $set->name(), 'sillyname';
ok $set->desc, 'a desc';
ok $set->accession, 'fakeaccesssion';

$res = Bio::Tools::HMMER::Results->new( -file => Bio::Root::IO->catfile("t","data","hmmsearch.out") , -type => 'hmmsearch');
my $seen =0;
ok $res->hmmfile, "HMM";
ok $res->seqfile, "HMM.dbtemp.29591";

my $first = 0;
foreach $set ( $res->each_Set) {
    foreach $domain ( $set->each_Domain ) {
    #print STDERR "Got domain ",$domain->seqname," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
      $seen = 1;
  }
}
ok $seen, 1;

ok $res->number, 1215, "\nBad number of domains. Expecting 1215. Got" . $res->number;

$res = Bio::Tools::HMMER::Results->new( -file => 
				      Bio::Root::IO->catfile("t","data",
							     "hmmpfam.out") , 
					-type => 'hmmpfam');

ok ($res->number, 2);

# parse HMM 2.2 files

$res = Bio::Tools::HMMER::Results->new( -file => 
				      Bio::Root::IO->catfile("t","data",
							     "L77119.hmmer"),
					-type => 'hmmpfam');
$seen =0;
ok $res->hmmfile, 'Pfam';
ok $res->seqfile, 'L77119.faa';
foreach $set ( $res->each_Set) {
    # only one set anyways

    ok($set->name, 'gi|1522636|gb|AAC37060.1|');
    ok($set->desc, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]');
    ok($set->accession, '[none]');
    foreach $domain ( $set->each_Domain ) {
	#print STDERR "Got domain ",$domain->seqname," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
	ok($domain->start, 280);
	ok($domain->end, 481);
	ok($domain->bits, -105.2);
	ok($domain->evalue, 0.0022 );
    }
}
ok ($res->number, 1);

# test for bugs #(1189,1034,1172)
$res = Bio::Tools::HMMER::Results->new( -file => Bio::Root::IO->catfile
					("t","data","hmmsearch.out") , 
					-type => 'hmmsearch');
my $res2 = $res->filter_on_cutoff(100,50);
ok($res2);
ok($res2->number, 604);


# now let's test the new Bio::SearchIO::hmmer
