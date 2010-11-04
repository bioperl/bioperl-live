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
    is(ref($result),'Bio::Search::Result::HMMERResult', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMPFAM', 'Check algorithm');
    is($result->algorithm_version, '2.1.1', 'Check algorithm version');
    is($result->hmm_name, 'pfam', 'Check hmm_name');
    is($result->sequence_file, '/home/birney/src/wise2/example/road.pep', 'Check sequence_file');
    is($result->query_name, 'roa1_drome', 'Check query_name');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 2, 'Check num_hits');
    my ($hsp,$hit);
    if( $hit = $result->next_model ) {
	is($hit->name, 'SEED', 'Check hit name');
	is($hit->raw_score, '146.1', 'Check hit raw_score');
	float_is($hit->significance, 6.3e-40, 'Check hit significance');
	is(ref($hit), 'Bio::Search::Hit::HMMERHit', 'Check for the correct hit reference type');
	is($hit->num_hsps, 1, 'Check num_hsps');

	if( defined( $hsp = $hit->next_domain ) ) {
	    is($hsp->hit->start, 1, 'Check for hit hmmfrom value');
	    is($hsp->hit->end, 77, 'Check for hit hmm to value');
	    is($hsp->query->start, 33, 'Check for query alifrom value');
	    is($hsp->query->end, 103, 'Check for query ali to value');
	    is($hsp->score, 71.2, 'Check for hsp score');
	    float_is($hsp->evalue, 2.2e-17, 'Check for hsp c-Evalue');
	    is($hsp->query_string, 'LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP',
                'Check for query string');
	    is($hsp->gaps('query'), 7, 'Check for number of gaps in query');
	    is($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG-kelggrklrv',
                'Check for hit string');
	    is($hsp->homology_string, 'lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ',
                'Check for homology string');
	    is(	length($hsp->homology_string), length($hsp->hit_string),
                'Check if homology string and hit string have an equal lenght');
	    is( length($hsp->query_string), length($hsp->homology_string),
                'Check if query string and homology string have an equal lenght');
	}
    }
    if( defined ($hit = $result->next_model) ) {
	if( defined($hsp = $hit->next_domain) ) {
	    is($hsp->hit->start, 1, 'Check for hit hmmfrom value');
	    is($hsp->hit->end, 77, 'Check for hit hmm to value');
	    is($hsp->query->start, 124, 'Check for query alifrom value');
	    is($hsp->query->end, 194, 'Check for query ali to value');
	    is($hsp->score, 75.5, 'Check for hsp score');
	    float_is($hsp->evalue, 1.1e-18, 'Check for hsp c-Evalue');
	    is($hsp->query_string, 'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV',
                'Check for query string');
	    is($hsp->gaps('query'), 6, 'Check for number of gaps in query');
	    is($hsp->hit_string, 'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv',
                'Check for hit string');
	    is($hsp->homology_string, 'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v',
                'Check for homology string');
	    is(	length($hsp->homology_string), length($hsp->hit_string),
                'Check if homology string and hit string have an equal lenght');
	    is( length($hsp->query_string), length($hsp->homology_string),
                'Check if query string and homology string have an equal lenght');
	}
	last;
    }
}
$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('hmmsearch.out'));
while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMSEARCH', 'Check algorithm');
    is($result->algorithm_version, '2.0', 'Check algorithm version');
    is($result->hmm_name, 'HMM [SEED]', 'Check hmm_name');
    is($result->sequence_file, 'HMM.dbtemp.29591', 'Check sequence_file');
    is($result->database_name, 'HMM.dbtemp.29591', 'Check database_name');
    is($result->query_name, 'SEED', 'Check query_name');
    is($result->query_description, '', 'Check query_description');
    is($result->num_hits(), 1215, 'Check num_hits');
    my $hit = $result->next_model;
    is($hit->name, 'Q91581', 'Check hit name');
    is($hit->description,'Q91581 POLYADENYLATION FACTOR 64 KDA SUBUN', 'Check for hit description');
    float_is($hit->significance, 2e-31, 'Check hit significance');
    is($hit->raw_score, 119.7, 'Check hit raw_score');
    my $hsp = $hit->next_domain;
    is($hsp->score,119.7, 'Check for hsp score');
    float_is($hsp->evalue, 2e-31, 'Check for hsp c-Evalue');
    is($hsp->query->start, 18, 'Check for query alifrom value');
    is($hsp->query->end, 89, 'Check for query ali to value');
    is($hsp->hit->start, 1, 'Check for hit hmmfrom value');
    is($hsp->hit->end, 77, 'Check for hit hmm to value');
    is($hsp->query->seq_id(), 'SEED', 'Check for query seq_id');
    is($hsp->hit->seq_id(), 'Q91581', 'Check for hit seq_id');
}

$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('L77119.hmmer'));

while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMPFAM', 'Check algorithm');
    is($result->algorithm_version, '2.2g', 'Check algorithm version');
    is($result->hmm_name, 'Pfam', 'Check hmm_name');
    is($result->sequence_file, 'L77119.faa', 'Check sequence_file');
    is($result->query_name, 'gi|1522636|gb|AAC37060.1|', 'Check query_name');
    is($result->query_description, 'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]',
        'Check query_description');
    is($result->num_hits(), 1, 'Check num_hits');
    my $hit = $result->next_hit;
    is($hit->name, 'Methylase_M', 'Check hit name');
    is($hit->description,'Type I restriction modification system, M', 'Check for hit description');
    float_is($hit->significance, 0.0022, 'Check hit significance');
    is($hit->raw_score, -105.2, 'Check hit raw_score');
    my $hsp = $hit->next_hsp;
    is($hsp->score,-105.2, 'Check for hsp score');
    float_is($hsp->evalue, 0.0022, 'Check for hsp evalue');
    is($hsp->query->start, 280, 'Check for query alifrom value');
    is($hsp->query->end, 481, 'Check for query ali to value');
    is($hsp->hit->start, 1, 'Check for hit hmmfrom value');
    is($hsp->hit->end, 279, 'Check for hit hmm to value');
    is($hsp->query->seq_id(), 'gi|1522636|gb|AAC37060.1|', 'Check for query seq_id');
    is($hsp->hit->seq_id(), 'Methylase_M', 'Check for hit seq_id');
    is($hsp->hit_string, 'lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG-dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn',
        'Check for hiy string');
    is($hsp->query_string, 'NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST',
        'Check for query string');
    is($hsp->homology_string, ' ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   ++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf ++++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+',
        'Check for homology string');
    is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '280 288 289 293-295 300 304 311 313-315 317 324-326 332 335 337 344-346 348 355 358-361 364-366 372 379 383-385 389 396 400 404-408 412 416 417 422 426 429-431 434-436 439 441 446 450 451 455 459 460 463 464 468 471 472 478',
        'Check for nomatch indices in query');
    is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '1 9 10 14-16 18-31 35 39 42-47 51-59 61 63-65 67 72-74 77-79 82 86 89-94 96 103-105 107 110 111 116 118 120-123 126-131 133 135-141 145 150 151 154 158-160 164 171 175 179-183 187 191-193 198 202 205-207 210-212 215 217 222 226 227 231 233 236 237 240-257 261 264-267 273 275-278',
        'Check for nomatch indices in hit');
	is(join(' ', $hsp->seq_inds('query', 'gap',1)), '296 306 309 321 328 334 335 350 356 366-368 376 417 456 463 470 479',
        'Check for gap indices in query');
	is(join(' ', $hsp->seq_inds('hit', 'gap',1)), '', 'Check for gap indices in hit');
}


$searchio = Bio::SearchIO->new(-format => 'hmmer',
			      -file   => test_input_file('cysprot1b.hmmsearch'));


while( my $result = $searchio->next_result ) {
    is(ref($result),'Bio::Search::Result::HMMERResult', 'Check for the correct result reference type');
    is($result->algorithm, 'HMMSEARCH', 'Check algorithm');
    is($result->algorithm_version, '2.2g', 'Check algorithm version');
    is($result->hmm_name, 'Peptidase_C1.hmm [Peptidase_C1]', 'Check hmm_name');
    is($result->database_name, 'cysprot1b.fa', 'Check database_name');
    is($result->sequence_file, 'cysprot1b.fa', 'Check sequence_file');
    is($result->query_name, 'Peptidase_C1', 'Check query_name');
    is($result->query_accession, 'PF00112', 'Check query_accession');
    is($result->query_description, 'Papain family cysteine protease', 'Check query_description');
    is($result->num_hits(), 4, 'Check num_hits');
    my $hit = $result->next_hit;
    is($hit->name, 'CATL_RAT', 'Check hit name');
    is($hit->description,'', 'Check for hit description');
    float_is($hit->significance, 2e-135, 'Check hit significance');
    is($hit->raw_score, 449.4, 'Check hit raw_score');
    my $hsp = $hit->next_hsp;
    is($hsp->score,449.4, 'Check for hsp score');
    float_is($hsp->evalue, 2e-135, 'Check for hsp evalue');
    is($hsp->query->start, 1, 'Check for query alifrom value');
    is($hsp->query->end, 337, 'Check for query ali to value');
    is($hsp->hit->start, 114, 'Check for hit hmmfrom value');
    is($hsp->hit->end, 332, 'Check for hit hmm to value');
    is($hsp->query->seq_id(), 'Peptidase_C1', 'Check for query seq_id');
    is($hsp->hit->seq_id(), 'CATL_RAT', 'Check for hit seq_id');
    is($hsp->hit_string, 'IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQMFLKT------GKLISLSEQNLVDCSH-DQGNQ------GCNG-GLMDFAFQYIKE-----NGGLDSEESY-----PYE----AKD-------------------GSCKYR-AEYAV-----ANDTGFVDIPQQ-----EKALMKAVATVGPISVAMDASHPS---LQFYSSG-------IYYEP---NCSSK---DLDHGVLVVGYGYEG-T------------------------------------DSNKDKYWLVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI',
        'Check for hiy string');
    is($hsp->homology_string, '+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+ ++kt      gkl+sLSEQ+LvDC++ d gn+      GCnG Glmd Af+Yik+     NgGl++E++Y     PY+    +kd                   g+Cky+  + ++     a+++g++d+p++     E+al+ka+a++GP+sVa+das+ s    q+Y+sG       +Y+++    C+++   +LdH+Vl+VGYG e+                                      ++++ +YW+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi',
        'Check for homology string');
    is($hsp->query_string, 'lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgtkawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikkeqIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgtCkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVaidasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGYGteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYWIVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi',
        'Check for query string');
    $hit = $result->next_hit;
    is($hit->name, 'CATL_HUMAN', 'Check hit name');
    is($hit->description,'', 'Check for hit description');
    float_is($hit->significance, 6.1e-134, 'Check hit significance');
    is($hit->raw_score, 444.5, 'Check hit raw_score');
}

# test for bug 2632 - CS lines should get ignored without breaking the parser
$searchio = Bio::SearchIO->new(-format => 'hmmer', -file => test_input_file('hmmpfam_cs.out'));
my $result = $searchio->next_result;
my $hit = $result->next_hit;
my $hsp = $hit->next_hsp;
is($hsp->seq_str, 'CGV-GFIADVNNVANHKIVVQALEALTCMEHRGACSADRDSGDGAGITTAIPWNLFQKSLQNQNIKFEQnDSVGVGMLFLPAHKLKES--KLIIETVLKEENLEIIGWRLVPTVQEVLGKQAYLNKPHVEQVFCKSSNLSKDRLEQQLFLVRKKIEKYIGINGKDwaheFYICSLSCYTIVYKGMMRSAVLGQFYQDLYHSEYTSSFAIYHRRFSTNTMPKWPLAQPMR---------FVSHNGEINTLLGNLNWMQSREPLLQSKVWKDRIHELKPITNKDNSDSANLDAAVELLIASGRSPEEALMILVPEAFQNQPDFA-NNTEISDFYEYYSGLQEPWDGPALVVFTNGKV-IGATLDRNGL-RPARYVIT----KDNLVIVSSES',
    'Check for hsp seq_str');
