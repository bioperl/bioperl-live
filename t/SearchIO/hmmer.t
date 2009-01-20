# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_hmmer.t 14989 2008-11-11 19:52:02Z cjfields $

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 116);
	
	use_ok('Bio::SearchIO');
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
	float_is($hit->significance, 6.3e-40);
	is(ref($hit), 'Bio::Search::Hit::HMMERHit');
	is($hit->num_hsps, 1);

	if( defined( $hsp = $hit->next_domain ) ) {
	    is($hsp->hit->start, 1);
	    is($hsp->hit->end, 77);
	    is($hsp->query->start, 33);
	    is($hsp->query->end, 103);
	    is($hsp->score, 71.2);
	    float_is($hsp->evalue, 2.2e-17);
	    is($hsp->query_string, 'LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP');
	    is($hsp->gaps('query'), 7);
	    is($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG-kelggrklrv');
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
	    float_is($hsp->evalue, 1.1e-18);
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
    float_is($hit->significance, 2e-31);
    is($hit->raw_score, 119.7);
    my $hsp = $hit->next_domain;
    is($hsp->score,119.7);
    float_is($hsp->evalue, 2e-31);
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
    float_is($hit->significance, 0.0022);
    is($hit->raw_score, -105.2);
    my $hsp = $hit->next_hsp;
    is($hsp->score,-105.2);
    float_is($hsp->evalue, 0.0022);
    is($hsp->query->start, 280);
    is($hsp->query->end, 481);
    is($hsp->hit->start, 1);
    is($hsp->hit->end, 279);
    is($hsp->query->seq_id(), 'gi|1522636|gb|AAC37060.1|');
    is($hsp->hit->seq_id(), 'Methylase_M');
    is($hsp->hit_string, 'lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG-dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn');
    is($hsp->query_string, 'NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST');
    is($hsp->homology_string, ' ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   ++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf ++++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+');
    is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '280 288 289 293-295 300 304 311 313-315 317 324-326 332 335 337 344-346 348 355 358-361 364-366 372 379 383-385 389 396 400 404-408 412 416 417 422 426 429-431 434-436 439 441 446 450 451 455 459 460 463 464 468 471 472 478');
    is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '1 9 10 14-16 18-31 35 39 42-47 51-59 61 63-65 67 72-74 77-79 82 86 89-94 96 103-105 107 110 111 116 118 120-123 126-131 133 135-141 145 150 151 154 158-160 164 171 175 179-183 187 191-193 198 202 205-207 210-212 215 217 222 226 227 231 233 236 237 240-257 261 264-267 273 275-278');
	is(join(' ', $hsp->seq_inds('query', 'gap',1)), '296 306 309 321 328 334 335 350 356 366-368 376 417 456 463 470 479');
	is(join(' ', $hsp->seq_inds('hit', 'gap',1)), '');
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
    float_is($hit->significance, 2e-135);
    is($hit->raw_score, 449.4);
    my $hsp = $hit->next_hsp;
    is($hsp->score,449.4);
    float_is($hsp->evalue, 2e-135);
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
    float_is($hit->significance, 6.1e-134);
    is($hit->raw_score, 444.5);
}

# test for bug 2632 - CS lines should get ignored without breaking the parser
$searchio = Bio::SearchIO->new(-format => 'hmmer', -file => test_input_file('hmmpfam_cs.out'));
my $result = $searchio->next_result;
my $hit = $result->next_hit;
my $hsp = $hit->next_hsp;
is $hsp->seq_str, 'CGV-GFIADVNNVANHKIVVQALEALTCMEHRGACSADRDSGDGAGITTAIPWNLFQKSLQNQNIKFEQnDSVGVGMLFLPAHKLKES--KLIIETVLKEENLEIIGWRLVPTVQEVLGKQAYLNKPHVEQVFCKSSNLSKDRLEQQLFLVRKKIEKYIGINGKDwaheFYICSLSCYTIVYKGMMRSAVLGQFYQDLYHSEYTSSFAIYHRRFSTNTMPKWPLAQPMR---------FVSHNGEINTLLGNLNWMQSREPLLQSKVWKDRIHELKPITNKDNSDSANLDAAVELLIASGRSPEEALMILVPEAFQNQPDFA-NNTEISDFYEYYSGLQEPWDGPALVVFTNGKV-IGATLDRNGL-RPARYVIT----KDNLVIVSSES';
