# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 140);
	
	use_ok('Bio::SearchIO');
	use_ok('Bio::Tools::HMMER::Domain');
	use_ok('Bio::Tools::HMMER::Set');
	use_ok('Bio::Tools::HMMER::Results');
}

my $searchio = Bio::SearchIO->new(-format => 'hmmer',
				 -file   => test_input_file('hmmpfam.out'));

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult');
    is($result->algorithm, 'HMMPFAM');
    is($result->algorithm_version, '2.1.1');
    is($result->hmm_name, 'pfam');
    is($result->sequence_file, '/home/birney/src/wise2/example/road.pep');
    is($result->query_name, 'roa1_drome');
    is($result->query_description, '');
    is($result->num_hits(), 2);
    my ($hsp,$hit);
    if( $hit = $result->next_model ) {
	is($hit->name, 'SEED');
	is($hit->raw_score, '146.1');
	is($hit->significance, '6.3e-40');
	is(ref($hit), 'Bio::Search::Hit::HMMERHit');
	is($hit->num_hsps, 1);

	if( defined( $hsp = $hit->next_domain ) ) {
	    is($hsp->hit->start, 1);
	    is($hsp->hit->end, 77);
	    is($hsp->query->start, 33);
	    is($hsp->query->end, 103);
	    is($hsp->score, 71.2);
	    is($hsp->evalue, '2.2e-17');
	    is($hsp->query_string, 'LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP');
	    is($hsp->gaps('query'), 7);
	    is($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG.kelggrklrv');
	    is($hsp->homology_string, 'lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ');
	    is(	length($hsp->homology_string), length($hsp->hit_string));
	    is( length($hsp->query_string), length($hsp->homology_string));
	}
    }
    if( defined ($hit = $result->next_model) ) {
	if( defined($hsp = $hit->next_domain) ) {
	    is($hsp->hit->start, 1);
	    is($hsp->hit->end, 77);
	    is($hsp->query->start, 124);
	    is($hsp->query->end, 194);
	    is($hsp->score, 75.5);
	    is($hsp->evalue, '1.1e-18');
	    is($hsp->query_string, 'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV');
	    is($hsp->gaps('query'), 6);
	    is($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv');	 
	    is($hsp->homology_string, 'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v');
	    is(	length($hsp->homology_string), length($hsp->hit_string));
	    is( length($hsp->query_string), length($hsp->homology_string));
	}
	last;
    }
}
$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('hmmsearch.out'));
while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult');
    is($result->algorithm, 'HMMSEARCH');
    is($result->algorithm_version, '2.0');
    is($result->hmm_name, 'HMM [SEED]');
    is($result->sequence_file, 'HMM.dbtemp.29591');
    is($result->database_name, 'HMM.dbtemp.29591');
    is($result->query_name, 'SEED');
    is($result->query_description, '');
    is($result->num_hits(), 1215);
    my $hit = $result->next_model;
    is($hit->name, 'Q91581');
    is($hit->description,'Q91581 POLYADENYLATION FACTOR 64 KDA SUBUN');
    is($hit->significance, '2e-31');
    is($hit->raw_score, 119.7);
    my $hsp = $hit->next_domain;
    is($hsp->score,119.7);
    is($hsp->evalue, '2e-31');
    is($hsp->query->start, 18);
    is($hsp->query->end, 89);
    is($hsp->hit->start, 1);
    is($hsp->hit->end, 77);
    is($hsp->query->seq_id(), 'SEED');
    is($hsp->hit->seq_id(), 'Q91581');   
}

$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('L77119.hmmer'));

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult');
    is($result->algorithm, 'HMMPFAM');
    is($result->algorithm_version, '2.2g');
    is($result->hmm_name, 'Pfam');
    is($result->sequence_file, 'L77119.faa');
    is($result->query_name, 'gi|1522636|gb|AAC37060.1|');
    is($result->query_description, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]');
    is($result->num_hits(), 1);
    my $hit = $result->next_hit;
    is($hit->name, 'Methylase_M');
    is($hit->description,'Type I restriction modification system, M');
    is($hit->significance, '0.0022');
    is($hit->raw_score, -105.2);
    my $hsp = $hit->next_hsp;
    is($hsp->score,-105.2);
    is($hsp->evalue, '0.0022');
    is($hsp->query->start, 280);
    is($hsp->query->end, 481);
    is($hsp->hit->start, 1);
    is($hsp->hit->end, 279);
    is($hsp->query->seq_id(), 'gi|1522636|gb|AAC37060.1|');
    is($hsp->hit->seq_id(), 'Methylase_M');
    is($hsp->hit_string, 'lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG.dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn');
    is($hsp->query_string, 'NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST');
    is($hsp->homology_string, ' ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   ++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf ++++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+');
    
}


$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('cysprot1b.hmmsearch'));


while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult');
    is($result->algorithm, 'HMMSEARCH');
    is($result->algorithm_version, '2.2g');
    is($result->hmm_name, 'Peptidase_C1.hmm [Peptidase_C1]');
    is($result->database_name, 'cysprot1b.fa');
    is($result->sequence_file, 'cysprot1b.fa');
    is($result->query_name, 'Peptidase_C1');
    is($result->query_accession, 'PF00112');
    is($result->query_description, 'Papain family cysteine protease');
    is($result->num_hits(), 4);
    my $hit = $result->next_hit;
    is($hit->name, 'CATL_RAT');
    is($hit->description,'');
    is($hit->significance, '2e-135');
    is($hit->raw_score, 449.4);
    my $hsp = $hit->next_hsp;
    is($hsp->score,449.4);
    is($hsp->evalue, '2e-135');
    is($hsp->query->start, 1);
    is($hsp->query->end, 337);
    is($hsp->hit->start, 114);
    is($hsp->hit->end, 332);
    is($hsp->query->seq_id(), 'Peptidase_C1');
    is($hsp->hit->seq_id(), 'CATL_RAT');
    is($hsp->hit_string, 'IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQMFLKT------GKLISLSEQNLVDCSH-DQGNQ------GCNG-GLMDFAFQYIKE-----NGGLDSEESY-----PYE----AKD-------------------GSCKYR-AEYAV-----ANDTGFVDIPQQ-----EKALMKAVATVGPISVAMDASHPS---LQFYSSG-------IYYEP---NCSSK---DLDHGVLVVGYGYEG-T------------------------------------DSNKDKYWLVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI');
    is($hsp->homology_string, '+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+ ++kt      gkl+sLSEQ+LvDC++ d gn+      GCnG Glmd Af+Yik+     NgGl++E++Y     PY+    +kd                   g+Cky+  + ++     a+++g++d+p++     E+al+ka+a++GP+sVa+das+ s    q+Y+sG       +Y+++    C+++   +LdH+Vl+VGYG e+                                      ++++ +YW+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi');
    is($hsp->query_string, 'lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgtkawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikkeqIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgtCkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVaidasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGYGteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYWIVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi');
    $hit = $result->next_hit;
    is($hit->name, 'CATL_HUMAN');
    is($hit->description,'');
    is($hit->significance, '6.1e-134');
    is($hit->raw_score, 444.5);
}
			      
my ($domain,$set,$homol,$rev,$res,$dom,@doms);
    $domain = Bio::Tools::HMMER::Domain->new(-verbose=>1);

is ref($domain), 'Bio::Tools::HMMER::Domain';

$domain->start(50);
$domain->end(200);
$domain->hstart(10);
$domain->hend(100);
$domain->seqbits(50);
$domain->bits(20);
$domain->evalue(0.0001);
$domain->seq_id('silly');


# test that we can get out forward and reverse homol_SeqFeatures
$homol = $domain->feature2();
is $homol->start(), 10;

$rev = $domain;

is $rev->start(), 50;

$set = Bio::Tools::HMMER::Set->new();
$set->add_Domain($domain);

@doms = $set->each_Domain();
$dom = shift @doms;

is $dom->start(), 50;

$set->bits(300);
$set->evalue(0.0001);
$set->name('sillyname');
$set->desc('a desc');
$set->accession('fakeaccesssion');
is $set->bits(), 300;
is $set->evalue(), 0.0001;
is $set->name(), 'sillyname';
is $set->desc, 'a desc';
is $set->accession, 'fakeaccesssion';

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmsearch.out') , -type => 'hmmsearch');
my $seen =0;
is $res->hmmfile, "HMM";
is $res->seqfile, "HMM.dbtemp.29591";

my $first = 0;
foreach $set ( $res->each_Set) {
    foreach $domain ( $set->each_Domain ) {
    #print STDERR "Got domain ",$domain->seq_id," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
      $seen = 1;
  }
}
is $seen, 1;

is $res->number, 1215;

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmpfam.out') , 
					-type => 'hmmpfam');

is ($res->number, 2);

# parse HMM 2.2 files

$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('L77119.hmmer'),
					-type => 'hmmpfam');
$seen =0;
is $res->hmmfile, 'Pfam';
is $res->seqfile, 'L77119.faa';
foreach $set ( $res->each_Set) {
    # only one set anyways

    is($set->name, 'gi|1522636|gb|AAC37060.1|');
    is($set->desc, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]');
    is($set->accession, '[none]');
    foreach $domain ( $set->each_Domain ) {
	#print STDERR "Got domain ",$domain->seq_id," start ",$domain->start," end ",$domain->end,"\n";
    # do nothing for the moment
	is($domain->start, 280);
	is($domain->end, 481);
	is($domain->bits, -105.2);
	is($domain->evalue, 0.0022 );
    }
}
is ($res->number, 1);

# test for bugs #(1189,1034,1172)
$res = Bio::Tools::HMMER::Results->new( -file => test_input_file('hmmsearch.out') , 
					-type => 'hmmsearch');
my $res2 = $res->filter_on_cutoff(100,50);
ok($res2);
is($res2->number, 604);

# now let's test the new Bio::SearchIO::hmmer (I believe this is in SearchIO.t)

