# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 327 );

    use_ok('Bio::SearchIO');
}

my $searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('hmmpfam.out')
);
my $result;

while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMPFAM', 'Check algorithm' );
    is( $result->algorithm_version, '2.1.1',   'Check algorithm version' );
    is( $result->hmm_name,          'pfam',    'Check hmm_name' );
    is( $result->sequence_file,
        '/home/birney/src/wise2/example/road.pep',
        'Check sequence_file'
    );
    is( $result->query_name,        'roa1_drome', 'Check query_name' );
    is( $result->query_description, '',           'Check query_description' );
    is( $result->num_hits(),        2,            'Check num_hits' );
    my ( $hsp, $hit );

    if ( $hit = $result->next_model ) {
        is( $hit->name,      'SEED',  'Check hit name' );
        is( $hit->raw_score, '146.1', 'Check hit raw_score' );
        float_is( $hit->significance, 6.3e-40, 'Check hit significance' );
        is( ref($hit), 'Bio::Search::Hit::HMMERHit',
            'Check for the correct hit reference type' );
        is( $hit->num_hsps, 1, 'Check num_hsps' );

        if ( defined( $hsp = $hit->next_domain ) ) {
            is( $hsp->hit->start,   1,    'Check for hit hmmfrom value' );
            is( $hsp->hit->end,     77,   'Check for hit hmm to value' );
            is( $hsp->query->start, 33,   'Check for query alifrom value' );
            is( $hsp->query->end,   103,  'Check for query ali to value' );
            is( $hsp->score,        71.2, 'Check for hsp score' );
            float_is( $hsp->evalue, 2.2e-17, 'Check for hsp c-Evalue' );
            is( $hsp->query_string,
                'LFIGGLDYRTTDENLKAHFEKWGNIVDVVVMKD-----PRTKRSRGFGFITYSHSSMIDEAQK--SRpHKIDGRVVEP',
                'Check for query string'
            );
            is( $hsp->gaps('query'), 7, 'Check for number of gaps in query' );
            is( $hsp->hit_string,
                'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnG-kelggrklrv',
                'Check for hit string'
            );
            is( $hsp->homology_string,
                'lf+g+L + +t+e Lk++F+k G iv++ +++D     + t++s+Gf+F+++  ++  + A +    +++++gr+++ ',
                'Check for homology string'
            );
            is( length( $hsp->homology_string ),
                length( $hsp->hit_string ),
                'Check if homology string and hit string have an equal lenght'
            );
            is( length( $hsp->query_string ),
                length( $hsp->homology_string ),
                'Check if query string and homology string have an equal lenght'
            );
        }
    }
    if ( defined( $hit = $result->next_model ) ) {
        if ( defined( $hsp = $hit->next_domain ) ) {
            is( $hsp->hit->start,   1,    'Check for hit hmmfrom value' );
            is( $hsp->hit->end,     77,   'Check for hit hmm to value' );
            is( $hsp->query->start, 124,  'Check for query alifrom value' );
            is( $hsp->query->end,   194,  'Check for query ali to value' );
            is( $hsp->score,        75.5, 'Check for hsp score' );
            float_is( $hsp->evalue, 1.1e-18, 'Check for hsp c-Evalue' );
            is( $hsp->query_string,
                'LFVGALKDDHDEQSIRDYFQHFGNIVDINIVID-----KETGKKRGFAFVEFDDYDPVDKVVL-QKQHQLNGKMVDV',
                'Check for query string'
            );
            is( $hsp->gaps('query'), 6, 'Check for number of gaps in query' );
            is( $hsp->hit_string,
                'lfVgNLppdvteedLkdlFskfGpivsikivrDiiekpketgkskGfaFVeFeseedAekAlealnGkelggrklrv',
                'Check for hit string'
            );
            is( $hsp->homology_string,
                'lfVg L  d +e+ ++d+F++fG iv+i+iv+D     ketgk +GfaFVeF++++ ++k +     ++l+g+ + v',
                'Check for homology string'
            );
            is( length( $hsp->homology_string ),
                length( $hsp->hit_string ),
                'Check if homology string and hit string have an equal lenght'
            );
            is( length( $hsp->query_string ),
                length( $hsp->homology_string ),
                'Check if query string and homology string have an equal lenght'
            );
        }
        last;
    }
}
$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('hmmsearch.out')
);
while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSEARCH',  'Check algorithm' );
    is( $result->algorithm_version, '2.0',        'Check algorithm version' );
    is( $result->hmm_name,          'HMM [SEED]', 'Check hmm_name' );
    is( $result->sequence_file, 'HMM.dbtemp.29591', 'Check sequence_file' );
    is( $result->database_name, 'HMM.dbtemp.29591', 'Check database_name' );
    is( $result->query_name,    'SEED',             'Check query_name' );
    is( $result->query_description, '',   'Check query_description' );
    is( $result->num_hits(),        1215, 'Check num_hits' );
    my $hit = $result->next_model;
    is( $hit->name, 'Q91581', 'Check hit name' );
    is( $hit->description,
        'Q91581 POLYADENYLATION FACTOR 64 KDA SUBUN',
        'Check for hit description'
    );
    float_is( $hit->significance, 2e-31, 'Check hit significance' );
    is( $hit->raw_score, 119.7, 'Check hit raw_score' );
    my $hsp = $hit->next_domain;
    is( $hsp->score, 119.7, 'Check for hsp score' );
    float_is( $hsp->evalue, 2e-31, 'Check for hsp c-Evalue' );
    is( $hsp->query->start,    18,       'Check for query alifrom value' );
    is( $hsp->query->end,      89,       'Check for query ali to value' );
    is( $hsp->hit->start,      1,        'Check for hit hmmfrom value' );
    is( $hsp->hit->end,        77,       'Check for hit hmm to value' );
    is( $hsp->query->seq_id(), 'SEED',   'Check for query seq_id' );
    is( $hsp->hit->seq_id(),   'Q91581', 'Check for hit seq_id' );
}

$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('L77119.hmmer')
);

while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMPFAM',    'Check algorithm' );
    is( $result->algorithm_version, '2.2g',       'Check algorithm version' );
    is( $result->hmm_name,          'Pfam',       'Check hmm_name' );
    is( $result->sequence_file,     'L77119.faa', 'Check sequence_file' );
    is( $result->query_name, 'gi|1522636|gb|AAC37060.1|',
        'Check query_name' );
    is( $result->query_description,
        'M. jannaschii predicted coding region MJECS02 [Methanococcus jannaschii]',
        'Check query_description'
    );
    is( $result->num_hits(), 1, 'Check num_hits' );
    my $hit = $result->next_hit;
    is( $hit->name, 'Methylase_M', 'Check hit name' );
    is( $hit->description,
        'Type I restriction modification system, M',
        'Check for hit description'
    );
    float_is( $hit->significance, 0.0022, 'Check hit significance' );
    is( $hit->raw_score, -105.2, 'Check hit raw_score' );
    my $hsp = $hit->next_hsp;
    is( $hsp->score, -105.2, 'Check for hsp score' );
    float_is( $hsp->evalue, 0.0022, 'Check for hsp evalue' );
    is( $hsp->query->start, 280, 'Check for query alifrom value' );
    is( $hsp->query->end,   481, 'Check for query ali to value' );
    is( $hsp->hit->start,   1,   'Check for hit hmmfrom value' );
    is( $hsp->hit->end,     279, 'Check for hit hmm to value' );
    is( $hsp->query->seq_id(),
        'gi|1522636|gb|AAC37060.1|', 'Check for query seq_id' );
    is( $hsp->hit->seq_id(), 'Methylase_M', 'Check for hit seq_id' );
    is( $hsp->hit_string,
        'lrnELentLWavADkLRGsmDaseYKdyVLGLlFlKYiSdkFlerrieieerktdtesepsldyakledqyeqlededlekedfyqkkGvFilPsqlFwdfikeaeknkldedigtdldkifseledqialgypaSeedfkGlfpdldfnsnkLgskaqarnetLtelidlfselelgtPmHNG-dfeelgikDlfGDaYEYLLgkFAeneGKsGGeFYTPqeVSkLiaeiLtigqpsegdfsIYDPAcGSGSLllqaskflgehdgkrnaisyYGQEsn',
        'Check for hiy string'
    );
    is( $hsp->query_string,
        'NTSELDKKKFAVLLMNR--------------LIFIKFLEDK------GIV---------PRDLLRRTYEDY---KKSNVLI-NYYDAY-L----KPLFYEVLNTPEDER--KENIRT-NPYYKDIPYL---N-G-------GLFRSNNV--PNELSFTIKDNEIIGEVINFLERYKFTLSTSEGsEEVELNP-DILGYVYEKLINILAEKGQKGLGAYYTPDEITSYIAKNT-IEPIVVE----------------RFKEIIK--NWKINDINF----ST',
        'Check for query string'
    );
    is( $hsp->homology_string,
        ' ++EL+++  av+   R              L+F K++ dk      +i+         p +   + +++y   ++   ++ ++y ++      + lF++++   e ++  ++++ + +    ++      + +       Glf ++++  ++ +s+   +ne ++e+i+ +++ +++     G++ +el   D++G +YE L+   Ae   K+ G +YTP e++  ia+ + i+  ++                  +++ ++    k+n+i +    s+',
        'Check for homology string'
    );
    is( join( ' ', $hsp->seq_inds( 'query', 'nomatch', 1 ) ),
        '280 288 289 293-295 300 304 311 313-315 317 324-326 332 335 337 344-346 348 355 358-361 364-366 372 379 383-385 389 396 400 404-408 412 416 417 422 426 429-431 434-436 439 441 446 450 451 455 459 460 463 464 468 471 472 478',
        'Check for nomatch indices in query'
    );
    is( join( ' ', $hsp->seq_inds( 'hit', 'nomatch', 1 ) ),
        '1 9 10 14-16 18-31 35 39 42-47 51-59 61 63-65 67 72-74 77-79 82 86 89-94 96 103-105 107 110 111 116 118 120-123 126-131 133 135-141 145 150 151 154 158-160 164 171 175 179-183 187 191-193 198 202 205-207 210-212 215 217 222 226 227 231 233 236 237 240-257 261 264-267 273 275-278',
        'Check for nomatch indices in hit'
    );
    is( join( ' ', $hsp->seq_inds( 'query', 'gap', 1 ) ),
        '296 306 309 321 328 334 335 350 356 366-368 376 417 456 463 470 479',
        'Check for gap indices in query'
    );
    is( join( ' ', $hsp->seq_inds( 'hit', 'gap', 1 ) ),
        '', 'Check for gap indices in hit' );
}

$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('cysprot1b.hmmsearch')
);

while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSEARCH', 'Check algorithm' );
    is( $result->algorithm_version, '2.2g',      'Check algorithm version' );
    is( $result->hmm_name,
        'Peptidase_C1.hmm [Peptidase_C1]',
        'Check hmm_name'
    );
    is( $result->database_name,   'cysprot1b.fa', 'Check database_name' );
    is( $result->sequence_file,   'cysprot1b.fa', 'Check sequence_file' );
    is( $result->query_name,      'Peptidase_C1', 'Check query_name' );
    is( $result->query_accession, 'PF00112',      'Check query_accession' );
    is( $result->query_description,
        'Papain family cysteine protease',
        'Check query_description'
    );
    is( $result->num_hits(), 4, 'Check num_hits' );
    my $hit = $result->next_hit;
    is( $hit->name,        'CATL_RAT', 'Check hit name' );
    is( $hit->description, '',         'Check for hit description' );
    float_is( $hit->significance, 2e-135, 'Check hit significance' );
    is( $hit->raw_score, 449.4, 'Check hit raw_score' );
    my $hsp = $hit->next_hsp;
    is( $hsp->score, 449.4, 'Check for hsp score' );
    float_is( $hsp->evalue, 2e-135, 'Check for hsp evalue' );
    is( $hsp->query->start, 1,   'Check for query alifrom value' );
    is( $hsp->query->end,   337, 'Check for query ali to value' );
    is( $hsp->hit->start,   114, 'Check for hit hmmfrom value' );
    is( $hsp->hit->end,     332, 'Check for hit hmm to value' );
    is( $hsp->query->seq_id(), 'Peptidase_C1', 'Check for query seq_id' );
    is( $hsp->hit->seq_id(),   'CATL_RAT',     'Check for hit seq_id' );
    is( $hsp->hit_string,
        'IPKTVDWRE-KG-CVTPVKNQG-QCGSCWAFSASGCLEGQMFLKT------GKLISLSEQNLVDCSH-DQGNQ------GCNG-GLMDFAFQYIKE-----NGGLDSEESY-----PYE----AKD-------------------GSCKYR-AEYAV-----ANDTGFVDIPQQ-----EKALMKAVATVGPISVAMDASHPS---LQFYSSG-------IYYEP---NCSSK---DLDHGVLVVGYGYEG-T------------------------------------DSNKDKYWLVKNSWGKEWGMDGYIKIAKDRN----NHCGLATAASYPI',
        'Check for hiy string'
    );
    is( $hsp->homology_string,
        '+P+++DWRe kg  VtpVK+QG qCGSCWAFSa g lEg+ ++kt      gkl+sLSEQ+LvDC++ d gn+      GCnG Glmd Af+Yik+     NgGl++E++Y     PY+    +kd                   g+Cky+  + ++     a+++g++d+p++     E+al+ka+a++GP+sVa+das+ s    q+Y+sG       +Y+++    C+++   +LdH+Vl+VGYG e+                                      ++++ +YW+VKNSWG++WG++GY++ia+++n    n+CG+a+ asypi',
        'Check for homology string'
    );
    is( $hsp->query_string,
        'lPesfDWReWkggaVtpVKdQGiqCGSCWAFSavgalEgryciktgtkawggklvsLSEQqLvDCdgedygnngesCGyGCnGGGlmdnAfeYikkeqIsnNgGlvtEsdYekgCkPYtdfPCgkdggndtyypCpgkaydpndTgtCkynckknskypktyakikgygdvpynvsTydEealqkalaknGPvsVaidasedskgDFqlYksGendvgyGvYkhtsageCggtpfteLdHAVliVGYGteneggtfdetssskksesgiqvssgsngssgSSgssgapiedkgkdYWIVKNSWGtdWGEnGYfriaRgknksgkneCGIaseasypi',
        'Check for query string'
    );
    $hit = $result->next_hit;
    is( $hit->name,        'CATL_HUMAN', 'Check hit name' );
    is( $hit->description, '',           'Check for hit description' );
    float_is( $hit->significance, 6.1e-134, 'Check hit significance' );
    is( $hit->raw_score, 444.5, 'Check hit raw_score' );
}

# test for bug 2632 - CS lines should get ignored without breaking the parser
$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('hmmpfam_cs.out')
);
$result = $searchio->next_result;
my $hit = $result->next_hit;
my $hsp = $hit->next_hsp;
is( $hsp->seq_str,
    'CGV-GFIADVNNVANHKIVVQALEALTCMEHRGACSADRDSGDGAGITTAIPWNLFQKSLQNQNIKFEQnDSVGVGMLFLPAHKLKES--KLIIETVLKEENLEIIGWRLVPTVQEVLGKQAYLNKPHVEQVFCKSSNLSKDRLEQQLFLVRKKIEKYIGINGKDwaheFYICSLSCYTIVYKGMMRSAVLGQFYQDLYHSEYTSSFAIYHRRFSTNTMPKWPLAQPMR---------FVSHNGEINTLLGNLNWMQSREPLLQSKVWKDRIHELKPITNKDNSDSANLDAAVELLIASGRSPEEALMILVPEAFQNQPDFA-NNTEISDFYEYYSGLQEPWDGPALVVFTNGKV-IGATLDRNGL-RPARYVIT----KDNLVIVSSES',
    'Check for hsp seq_str'
);

# Tests for hmmer3 output here
$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmscan.out'),
    -verbose => 1
);
is( ref($searchio), 'Bio::SearchIO::hmmer3',
    'Check if correct searchio object is returned' );
while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSCAN', 'Check algorithm' );
    is( $result->algorithm_version, '3.0',     'Check algorithm version' );
    is( $result->hmm_name,
        '/data/biodata/HMMerDB/Pfam.hmm',
        'Check hmm_name'
    );
    is( $result->sequence_file, 'BA000019.orf1.fasta',
        'Check sequence_file' );
    is( $result->query_name, 'BA000019.orf1', 'Check query_name' );
    is( $result->query_length, '198', 'Check query_length' );
    is( $result->query_description, '', 'Check query_description' );
    is( $result->num_hits(),        1,  'Check num_hits' );
    my ( $hsp, $hit );

    if ( $hit = $result->next_model ) {
        is( ref($hit), 'Bio::Search::Hit::HMMERHit',
            'Check for the correct hit reference type' );
        is( $hit->name, 'Peripla_BP_2', 'Check hit name' );
        is( $hit->description,
            'Periplasmic binding protein',
            'Check for hit description'
        );
        is( $hit->raw_score, '105.2', 'Check hit raw_score' );
        float_is( $hit->significance, 6e-30, 'Check hit significance' );
        is( $hit->num_hsps, 1, 'Check num_hsps' );

        if ( defined( $hsp = $hit->next_domain ) ) {
            is( ref($hsp), 'Bio::Search::HSP::HMMERHSP',
                'Check for correct hsp reference type' );
            is( $hsp->hit->start,   59,  'Check for hit hmmfrom value' );
            is( $hsp->hit->end,     236, 'Check for hit hmm to value' );
            is( $hsp->query->start, 2,   'Check for query alifrom value' );
            is( $hsp->query->end,   173, 'Check for query ali to value' );
            is( $hsp->score, '105.0', 'Check for hsp score' );
            float_is( $hsp->evalue, 1.5e-33, 'Check for hsp c-Evalue' );
            is( $hsp->query_string,
                'LKPDLIIGREYQ---KNIYNQLSNFAPTVLVDWGSF-TSFQDNFRYIAQVLNEEEQGKLVLQQYQKRIRDLQDRMGERlQKIEVSVIGFSGQSIKSLNR-DAVFNQVLDDAGIKRIsIQKNQQERYLEISIENLNKYDADVLFVINE---SKEQLYPDLKNPLWHHLRAVKKQQVYVVNQ',
                'Check for query string'
            );
            is( $hsp->hit_string,
                'lkPDlvivsafgalvseieellelgipvvavessstaeslleqirllgellgeedeaeelvaelesridavkaridsl-kpktvlvfgyadegikvvfgsgswvgdlldaaggeni-iaeakgseseeisaEqilaadpdviivsgrgedtktgveelkenplwaelpAvkngrvyllds',
                'Check for hit string'
            );
            is( $hsp->homology_string,
                'lkPDl+i+ +++   ++i+++l++ +p+v v+  s+  s+++ +r ++++l+ee++++ + +++++ri+++++r  +  ++ +v+v+g+++ +ik+++  +  ++++ld+ag++ i i++++++ + eis+E+++++d+dv++v       k+ +   ++nplw +l+Avk+++vy++++',
                'Check for homology string'
            );
        }
    }
}

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmsearch3.out'),
    -verbose => 1
);
while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSEARCH', 'Check algorithm' );
    is( $result->algorithm_version, '3.0',       'Check algorithm version' );
    is( $result->hmm_name,          'Kv9.hmm',   'Check hmm_name' );
    is( $result->sequence_file,
        '/home/pboutet/Desktop/databases/nr_May26',
        'Check sequence_file'
    );
    is( $result->query_name,        'Kv9', 'Check query_name' );
    is( $result->query_length,      '481', 'Check query_length' );
    is( $result->query_description, '',    'Check query_description' );
    is( $result->num_hits(),        2,     'Check num_hits' );

    while ( my $hit = $result->next_model ) {
    }
}

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmscan_multi_domain.out'),
    -verbose => 1
);

my @multi_hits = (
    [   'PPC',
        'Bacterial pre-peptidase C-terminal domain',
        '111.0', 3.1e-32, 6,
        [   [ 4,  59, 117,  183,  0.5,  0.16 ],
            [ 12, 58, 347,  388,  -0.6, 0.36 ],
            [ 1,  69, 470,  549,  71.3, 1.3e-23 ],
            [ 15, 25, 582,  603,  -3.2, 2 ],
            [ 13, 36, 987,  1019, -1.1, 0.5 ],
            [ 1,  69, 1087, 1168, 54.4, 2.4e-18 ]
        ]
    ],
    [   'HemolysinCabind',
        'Hemolysin-type calcium-binding repeat (2 cop',
        '47.9', 4.7e-13, 3,
        [   [ 2, 13, 1214, 1225, 5.9,  0.0026 ],
            [ 1, 18, 1231, 1248, 10.8, 6.8e-5 ],
            [ 4, 18, 1243, 1257, 11.4, 4.3e-05 ]
        ]
    ]
);

while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSCAN', 'Check algorithm' );
    is( $result->algorithm_version, '3.0',     'Check algorithm version' );
    is( $result->hmm_name,
        '/data/biodata/HMMerDB/Pfam-A.hmm',
        'Check hmm_name'
    );
    is( $result->sequence_file, 'BA000019.orf37.fasta',
        'Check sequence_file' );
    is( $result->query_name, 'BA000019.orf37', 'Check query_name' );
    is( $result->query_length, '1418', 'Check query_length' );
    is( $result->query_description, '', 'Check query_description' );
    is( $result->num_hits(),        2,  'Check num_hits' );
    my ( $hsp, $hit );

    while ( $hit = $result->next_model ) {
        my @expected = @{ shift @multi_hits };
        is( ref($hit), 'Bio::Search::Hit::HMMERHit',
            'Check for the correct hit reference type' );
        is( $hit->name,        shift @expected, 'Check hit name' );
        is( $hit->description, shift @expected, 'Check for hit description' );
        is( $hit->raw_score,   shift @expected, 'Check hit raw_score' );
        float_is(
            $hit->significance,
            shift @expected,
            'Check hit significance'
        );
        is( $hit->num_hsps, shift @expected, 'Check num_hsps' );
        my @hsp_list = @{ shift @expected };

        while ( defined( $hsp = $hit->next_domain ) ) {
            my @hsp_exp = @{ shift @hsp_list };
            is( ref($hsp), 'Bio::Search::HSP::HMMERHSP',
                'Check for correct hsp reference type' );
            is( $hsp->hit->start,
                shift @hsp_exp,
                'Check for hit envfrom value'
            );
            is( $hsp->hit->end, shift @hsp_exp,
                'Check for hit env to value' );
            is( $hsp->query->start,
                shift @hsp_exp,
                'Check for query hmmfrom value'
            );
            is( $hsp->query->end,
                shift @hsp_exp,
                'Check for query hmm to value'
            );
            is( $hsp->score, shift @hsp_exp, 'Check for hsp score' );
            float_is( $hsp->evalue, shift @hsp_exp,
                'Check for hsp c-Evalue' );
        }
    }
}

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmscan_sec_struct.out'),
    -verbose => 1
);

@multi_hits = (
    [   'HTH_AraC',
        'Bacterial regulatory helix-turn-helix proteins, Ara',
        '41.3', 6.7e-11, 2,
        [   [ 'siadiAeevgfSpsyfsrlFkkytGvt', 'SLMELSRQVGLNDCTLKRGFRLVFDTT' ],
            [   'nwsiadiAeevgf-SpsyfsrlFkkytGvtPsqyr',
                'EINISQAARRVGFsSRSYFATAFRKKFGINPKEFL'
            ]
        ]
    ],
    [   'PKSI-KS_m3',
        '', '38.2', 3.8e-12, 2,
        [   [ 'GPSvtVDTACSSSLvA', 'GPSVTVDTLCSSSLVA' ],
            [ 'GPSvtVDTACSSSLv',  'GPNLVIDSACSSALV' ]
        ]
    ],
    [   'DUF746',
        'Domain of Unknown Function (DUF746)',
        '13.9', 0.023, 2,
        [   [   'rllIrlLsqplslaeaadqlgtdegiiak',
                'EILIRNLENPPSLMELSRQVGLNDCTLKR'
            ],
            [ 'plslaeaadqlgtdeg', 'EINISQAARRVGFSSR' ]
        ]
    ]
);

while ( $result = $searchio->next_result ) {
    is( ref($result),
        'Bio::Search::Result::HMMERResult',
        'Check for the correct result reference type'
    );
    is( $result->algorithm,         'HMMSCAN',    'Check algorithm' );
    is( $result->algorithm_version, '3.0',        'Check algorithm version' );
    is( $result->hmm_name,          'Pfam-A.hmm', 'Check hmm_name' );
    is( $result->sequence_file, 'BA000019.orf8.fasta',
        'Check sequence_file' );
    is( $result->query_name, 'BA000019.orf8', 'Check query_name' );
    is( $result->query_length, '348', 'Check query_length' );
    is( $result->query_description, '', 'Check query_description' );
    is( $result->num_hits(),        3,  'Check num_hits' );
    my ( $hsp, $hit );

    while ( $hit = $result->next_model ) {
        my @expected = @{ shift @multi_hits };
        is( ref($hit), 'Bio::Search::Hit::HMMERHit',
            'Check for the correct hit reference type' );
        is( $hit->name,        shift @expected, 'Check hit name' );
        is( $hit->description, shift @expected, 'Check for hit description' );
        is( $hit->raw_score,   shift @expected, 'Check hit raw_score' );
        float_is(
            $hit->significance,
            shift @expected,
            'Check hit significance'
        );
        is( $hit->num_hsps, shift @expected, 'Check num_hsps' );
        my @hsp_list = @{ shift @expected };

        while ( defined( $hsp = $hit->next_domain ) ) {
            my @hsp_exp = @{ shift @hsp_list };
            is( ref($hsp), 'Bio::Search::HSP::HMMERHSP',
                'Check for correct hsp reference type' );
            is( $hsp->hit_string,   shift @hsp_exp, 'Check hit sequence' );
            is( $hsp->query_string, shift @hsp_exp, 'Check query sequence' );
        }
    }
}

# Make sure that you can also directly call the hmmer2 and hmmer3 subclasses
$searchio = Bio::SearchIO->new(
    -format => 'hmmer2',
    -file   => test_input_file('hmmpfam.out')
);
is( ref($searchio), 'Bio::SearchIO::hmmer2',
    'Check if loading hmmpfam output via the hmm2 parser directly works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format => 'hmmer2',
    -file   => test_input_file('hmmsearch.out')
);
is( ref($searchio), 'Bio::SearchIO::hmmer2',
    'Check if loading hmmsearch2 output via the hmm2 parser directly works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format => 'hmmer3',
    -file   => test_input_file('hmmscan.out')
);
is( ref($searchio), 'Bio::SearchIO::hmmer3',
    'Check if loading hmmscan output via the hmm3 parser directly works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -file   => test_input_file('hmmsearch3.out')
);
is( ref($searchio), 'Bio::SearchIO::hmmer3',
    'Check if loading hmmsearch3 output via the hmm3 parser directly works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

# Make sure that you can also specify the -version parameter directly
$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmpfam.out'),
    -version => 2
);
is( ref($searchio), 'Bio::SearchIO::hmmer2',
    'Check if selecting the correct hmmpfam parser using -version works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmsearch.out'),
    -version => 2
);
is( ref($searchio), 'Bio::SearchIO::hmmer2',
    'Check if selecting the correct hmmsearch2 parser using -version works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer3',
    -file    => test_input_file('hmmscan.out'),
    -version => 3
);
is( ref($searchio), 'Bio::SearchIO::hmmer3',
    'Check if selecting the correct hmmscan parser using -version works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

$searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -file    => test_input_file('hmmsearch3.out'),
    -version => 3
);
is( ref($searchio), 'Bio::SearchIO::hmmer3',
    'Check if selecting the correct hmmsearch3 parser using -version works' );
is( ref( $searchio->next_result ),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);

my $pipestr = "cat " . test_input_file('hmmpfam.out') . " |";
open( my $pipefh, $pipestr );

$searchio = Bio::SearchIO->new(
    -format => 'hmmer',
    -fh     => $pipefh
);
is( ref($searchio), 'Bio::SearchIO::hmmer2',
    'Check if reading from a pipe works' );
$result = $searchio->next_result;
is( ref($result),
    'Bio::Search::Result::HMMERResult',
    'Check for the correct result reference type'
);
is( $result->num_hits(), 2, 'Check num_hits' );

# bug 3376
{
    my $in = Bio::SearchIO->new(
        -format => 'hmmer',
        -file   => test_input_file('pfamOutput-bug3376.out')
    );
    my $result = $in->next_result;
    my $hit    = $result->next_hit;
    my $hsp    = $hit->next_hsp;
    is( $hsp->hit_string,
        "svfqqqqssksttgstvtAiAiAigYRYRYRAvtWnsGsLssGvnDnDnDqqsdgLYtiYYsvtvpssslpsqtviHHHaHkasstkiiikiePr",
        "bug3376"
    );
}
# end bug 3376

# bug 3421 - making sure a full line of dashes in an HSP is parsed correctly
{
    my $in = Bio::SearchIO->new(
        -format => 'hmmer',
        -file   => test_input_file('hmmpfam_HSPdashline.txt')
    );
    my $result = $in->next_result;
    my $hit    = $result->next_hit;
    my $hsp    = $hit->next_hsp;
    is( $hsp->length, "561",
        "bug3421 - Check if can correctly parse an HSP with line full of dashes"
    );
}
# end bug 3421

# bug 3302
{
    my $in = Bio::SearchIO->new(
        -format => 'hmmer',
        -file   => test_input_file('hmmpfam_multiresult.out')
    );
    my $result = $in->next_result;
    $result = $in->next_result;
    my $hit = $result->next_hit;
    is( $hit->name, "IS66_ORF3.uniq", "bug3302 - Check if can parse multiresult hmmer" );
}
# end bug 3302

# HMMER 3.1 nhmmer output
{
    my $in = Bio::SearchIO->new(
        -format  => 'hmmer',
        -version => 3,
        -file    => test_input_file('nhmmer-3.1.out')
    );
    my $result = $in->next_result;
    is( $result->algorithm_version, '3.1b1', 'Check nhmmer algorithm version' );

    my $hit = $result->next_hit;
    is( $hit->name, "seq1", "Check nhmmer hit name" );
    is( $hit->description, "Description of seq1", "Check nhmmer hit description" );
    is( $hit->significance, "3.2e-48", "Check nhmmer hit significance" );
    is( $hit->score,        "148.2",   "Check nhmmer hit score" );

    my $hsp = $hit->next_hsp;
    is( $hsp->score, "148.2",   "Check nhmmer hsp score" );
    is( $hsp->significance, "3.2e-48",   "Check nhmmer hsp score" );
    is( $hsp->start('hit'), "1",   "Check nhmmer hsp hit start" );
    is( $hsp->end('hit'), "151",   "Check nhmmer hsp hit end" );
    is( $hsp->start('query'), "258",   "Check nhmmer hsp query start" );
    is( $hsp->end('query'), "411",   "Check nhmmer hsp query end" );
    is( $hsp->strand('hit'), '1',   "Check nhmmer hsp hit strand" );
    is( $hsp->strand('query'), '1',   "Check nhmmer hsp query strand" );

    $hit = $result->next_hit;
    is( $hit->name, "seq2", "Check nhmmer hit name" );
    is( $hit->description, "Description of seq2", "Check nhmmer hit description" );
    is( $hit->significance, "3.9e-15", "Check nhmmer hit significance" );
    is( $hit->score,        "38.6",   "Check nhmmer hit score" );

    $hsp = $hit->next_hsp;
    is( $hsp->score, '38.6',   "Check nhmmer hsp score" );
    is( $hsp->significance, '3.9e-15',   "Check nhmmer hsp score" );
    is( $hsp->start('query'), '34',   "Check nhmmer hsp query start" );
    is( $hsp->end('query'), '92',   "Check nhmmer hsp query end" );
    is( $hsp->start('hit'), '1',   "Check nhmmer hsp hit start" );
    is( $hsp->end('hit'), '59',   "Check nhmmer hsp hit end" );
    is( $hsp->strand('hit'), '-1',   "Check nhmmer hsp hit strand" );
    is( $hsp->strand('query'), '1',   "Check nhmmer hsp query strand" );
}
# end HMMER 3.1 nhmmer output

