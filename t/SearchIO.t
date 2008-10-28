# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 1811);
	
	use_ok('Bio::SearchIO');
	use_ok('Bio::SearchIO::Writer::HitTableWriter');
	use_ok('Bio::SearchIO::Writer::HTMLResultWriter');
}

my ($searchio, $result,$iter,$hit,$hsp);

SKIP: {
	# XML encoding/decoding done within XML::SAX now, though some parsers
	# do not work properly (XML::SAX::PurePerl, XML::LibXML::SAX)
	test_skip(-tests => 129, -requires_module => 'XML::SAX');
	
    eval {
		# test with RPSBLAST data first
		# this needs to be eval'd b/c the XML::SAX parser object is
		# instantiated in the constructor
		$searchio = Bio::SearchIO->new('-tempfile' => 1,
			   '-format' => 'blastxml',
			   '-file'   => test_input_file('ecoli_domains.rps.xml'),
               '-blasttype' => 'blast',
			   '-verbose' => -1);
		# PurePerl works with these BLAST reports, so removed verbose promotion
		$result = $searchio->next_result;
		die if !defined $result;
	};
	if ($@ && $@ =~ m{Handler could not resolve external entity}) {
		skip("XML::SAX::Expat does not work with XML tests; skipping",129);
	} elsif ($@) {
		skip("Problem with XML::SAX setup: $@. Check ParserDetails.ini; skipping XML tests",129);
	}
	
    isa_ok($result, 'Bio::Search::Result::ResultI');
    is($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam', 'database_name()');
    is($result->query_name,'gi|1786182|gb|AAC73112.1|','query_name()');
    is($result->query_description, '(AE000111) thr operon leader peptide [Escherichia coli]');
    is($result->query_accession, 'AAC73112.1');
    is($result->query_gi, 1786182);
    is($result->query_length, 21);
    is($result->algorithm, 'BLASTP');
    is($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');

    is($result->available_parameters, 8);
    is($result->get_parameter('gapext'), 1);
    is($result->available_statistics, 5);
    is($result->get_statistic('lambda'), 0.267);

	# this result actually has a hit
    $result = $searchio->next_result;
    $hit = $result->next_hit;
    is($hit->name, 'gnl|Pfam|pfam00742');
    is($hit->description(), 'HomoS_dh, HomoS dehydrogenase');
    is($hit->accession, 'pfam00742');
    is($hit->ncbi_gi, ''); # not found
    is($hit->length, 310);

    $hsp = $hit->next_hsp;
    is($hsp->query->seq_id, $result->query_name,'query name on HSP');
    is($hsp->query->seqdesc, $result->query_description,'query desc on HSP');
    is($hsp->hit->seq_id, $hit->name,'hitname');
    is($hsp->hit->seqdesc, $hit->description,'hitdesc');
    is($hsp->pvalue, undef);
    is(sprintf("%g",$hsp->evalue), sprintf("%g",'1.46134e-90'));
    is($hsp->score, 838);
    is($hsp->bits,327.405);
    is($hsp->query->start, 498);
    is($hsp->query->end,815);
    is($hsp->hit->start, 3);
    is($hsp->hit->end, 310);
    is($hsp->query->frame,0);
    is($hsp->hit->frame,0);
    is(sprintf("%.2f", $hsp->percent_identity), 37.73);
    is(sprintf("%.4f", $hsp->frac_identical('hit')), 0.3994);
    is(sprintf("%.4f", $hsp->frac_identical('query')), 0.3868);
    is(sprintf("%.4f",$hsp->query->frac_identical), 0.3868);

    is(sprintf("%.4f",$hsp->frac_conserved('total')),0.5245);
    is(sprintf("%.4f",$hsp->frac_conserved('hit')),0.5552);
    is(sprintf("%.4f",$hsp->frac_conserved('query')),0.5377);
    is($hsp->gaps('total'), 26);
    is($hsp->length('hsp'), 326);
    is($hsp->query_string, 'LRVCGVANSKALLTNVHGLNLENWQEELAQAKEPF-NLGRLIRLVKEYHLLN----PVIVDCTSSQAVAD-QYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE-GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARET-GRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS');
    is($hsp->hit_string, 'GVVTGITDSREMLLSRIGLPLEIWKVALRDLEKPRKDLGKLDLTDDAFAVVDDPDIDVVVELTGGIEVARELYLDALEEGKHVVTANKALNASHGDEYLAL---AEKSGVDVLYEAAVAGGIPIIKTLRELLATGDRILKIEGIFNGTTNFILSEMDEKGLPFSDVLAEAQELGYTEADPRDDVEGIDAARKLAILARIAFGIELELDDVYVEGISPITAEDISSADEFGYTLKLLDEAMRQRVEDAESGGEVLRYPTLIPE-------------DHPLASVKGSDNAVAVEGEAYG--PLMFYGPGAGAEPTASAVVADIVRIAR');
    is($hsp->homology_string, '  V G+ +S+ +L +  GL LE W+  L   ++P  +LG+L      + +++     V+V+ T    VA   Y D L EG HVVT NK  N S  D Y  L   AEKS    LY+  V  G+P+I+ L+ LL  GD ++K  GI +G+ ++I  ++DE G+ FS+    A+E+GYTE DPRDD+ G+D ARKL ILAR   G ELEL D+ +E + P           F   L  LD+    RV  A   G+VLRY   I E             + PL  VK  +NA+A     Y   PL+  G GAG + TA+ V AD++R   ');
    is(join(' ', $hsp->seq_inds('query', 'gap',1)), '532 548 562 649 690');
    is($hsp->ambiguous_seq_inds, '');
	
	# one more 
    $hit = $result->next_hit;
    isa_ok($hit,'Bio::Search::Hit::HitI');
    
    my $results_left = 8;
    while( $result = $searchio->next_result ) { ok($result); $results_left--; }
    is($results_left, 0);


    $searchio = Bio::SearchIO->new(-format => 'blastxml',
								  -verbose => -1,
				  -file => test_input_file('plague_yeast.bls.xml'));

    $result = $searchio->next_result;

    is($result->database_name, 'yeast.aa');
    is($result->query_name, 'gi|5763811|emb|CAB53164.1|');
    is($result->query_description,  'putative transposase [Yersinia pestis]');
    is($result->query_accession, 'CAB53164.1');
    is($result->query_gi, 5763811);  
    is($result->query_length, 340);

    $hit = $result->next_hit;
    ok(! $hit);

    $searchio = Bio::SearchIO->new(-format => 'blastxml',
								  -verbose => -1,
				  -file => test_input_file('mus.bls.xml'));

    $result = $searchio->next_result;

    is($result->database_name,'Hs15_up1000');
    is($result->query_name,'NM_011441_up_1000_chr1_4505586_r');
    is($result->query_description,'chr1:4505586-4506585');
    is($result->query_accession,'NM_011441_up_1000_chr1_4505586_r');
    is($result->query_gi, '');
    is($result->query_length,'1000');
    $hit = $result->next_hit;
    is($hit->name,'NM_001938_up_1000_chr1_93161154_f');
    is($hit->description,'chr1:93161154-93162153');
    is($hit->ncbi_gi, ''); # none reported
    is($hit->accession,'3153');
    is($hit->length,'1000');
    
    # deal with new BLAST XML changes
    $searchio = Bio::SearchIO->new(-format => 'blastxml',
								  -verbose => -1,
				  -file => test_input_file('newblast.xml'));

    $result = $searchio->next_result;

    is($result->database_name,'nr');
    is($result->database_name,'nr');
    is($result->database_letters,'1479795817');
    is($result->database_entries,'4299737');
    is($result->algorithm,'BLASTP');
    is($result->algorithm_version,'BLASTP 2.2.15 [Oct-15-2006]');
    
    # some XML::SAX parsers (PurePerl, XML::SAX::LibXML) don't decode entities
    # properly, not fixable using decode_entities()
    like($result->algorithm_reference, qr{Nucleic Acids Res} ); 
    is($result->available_parameters,4);
    is($result->available_statistics,5);
    is($result->query_name,'gi|15600734|ref|NP_254228.1|');
    is($result->query_description,'dihydroorotase [Pseudomonas aeruginosa PAO1]');
    is($result->query_accession,'NP_254228.1');
    is($result->query_gi, 15600734);
    is($result->query_length,'445');
    $hit = $result->next_hit;
    is($hit->name,'gi|15600734|ref|NP_254228.1|');
    is($hit->description,'dihydroorotase [Pseudomonas aeruginosa PAO1] '.
       '>gi|107104643|ref|ZP_01368561.1| hypothetical protein PaerPA_01005722 '.
       '[Pseudomonas aeruginosa PACS2] >gi|9951880|gb|AAG08926.1|AE004966_8 '.
       'dihydroorotase [Pseudomonas aeruginosa PAO1]');
    is($hit->accession,'NP_254228');
    is($hit->length,'445');
    $hsp = $hit->next_hsp;
    is($hsp->query->seq_id, $result->query_name,'query name on HSP');
    is($hsp->query->seqdesc, $result->query_description,'query desc on HSP');
    is($hsp->hit->seq_id, $hit->name,'hitname');
    is($hsp->hit->seqdesc, $hit->description,'hitdesc');
    is($hsp->pvalue, undef);
    is(sprintf("%g",$hsp->evalue), sprintf("%g",'0'));
    is($hsp->score, 2251);
    is($hsp->bits,871.692);
    is($hsp->query->start, 1);
    is($hsp->query->end,445);
    is($hsp->hit->start, 1);
    is($hsp->hit->end, 445);
    is($hsp->query->frame,0);
    is($hsp->hit->frame,0);
    
    $result = $searchio->next_result;

    is($result->database_name,'nr'); 
    is($result->database_letters,'1479795817'); 
    is($result->database_entries,'4299737');
    is($result->algorithm,'BLASTP');
    is($result->algorithm_version,'BLASTP 2.2.15 [Oct-15-2006]'); 
    like($result->algorithm_reference, qr{Nucleic Acids Res} );
    is($result->available_parameters,4); 
    is($result->available_statistics,5);
    is($result->query_name,'gi|15598723|ref|NP_252217.1|');
    is($result->query_description,'dihydroorotase [Pseudomonas aeruginosa PAO1]');
    is($result->query_accession,'NP_252217.1');
    is($result->query_gi, 15598723);
    is($result->query_length,'348');
    $hit = $result->next_hit;
    is($hit->name,'gi|15598723|ref|NP_252217.1|');
    is($hit->description,'dihydroorotase [Pseudomonas aeruginosa PAO1] '.
       '>gi|6226683|sp|P72170|PYRC_PSEAE Dihydroorotase (DHOase) '.
       '>gi|9949676|gb|AAG06915.1|AE004773_4 dihydroorotase [Pseudomonas aeruginosa PAO1] '.
       '>gi|3868712|gb|AAC73109.1| dihydroorotase [Pseudomonas aeruginosa]');
    is($hit->ncbi_gi, 15598723);
    is($hit->accession,'NP_252217');
    is($hit->length,'348');
    $hsp = $hit->next_hsp;
    is($hsp->query->seq_id, $result->query_name,'query name on HSP');
    is($hsp->query->seqdesc, $result->query_description,'query desc on HSP');
    is($hsp->hit->seq_id, $hit->name,'hitname');
    is($hsp->hit->seqdesc, $hit->description,'hitdesc');
    is($hsp->pvalue, undef);
    is(sprintf("%g",$hsp->evalue), sprintf("%g",'0'));
    is($hsp->score, 1780);
    is($hsp->bits,690.263);
    is($hsp->query->start, 1);
    is($hsp->query->end,348);
    is($hsp->hit->start, 1);
    is($hsp->hit->end, 348);
    is($hsp->query->frame,0);
    is($hsp->hit->frame,0);
    
    # PSIBLAST XML parsing 
    
    $searchio = Bio::SearchIO->new('-tempfile' => 1,
           '-format' => 'blastxml',
           '-file'   => test_input_file('psiblast.xml'),
           '-blasttype' => 'psiblast');
    
    my $result = $searchio->next_result;
    is($result->database_name, 'AL591824.faa');
    is($result->database_entries, 2846);
    is($result->database_letters, 870878);
    is($result->algorithm, 'BLASTP');
    like($result->algorithm_version, qr/2\.2\.16/);
    is($result->query_name, 'gi|1373160|gb|AAB57770.1|');
    is($result->query_accession, 'AAB57770.1');
    is($result->query_gi, '1373160');
    is($result->query_length, 173);
    is($result->get_statistic('kappa') , 0.0475563);
    cmp_ok($result->get_statistic('lambda'), '==', 0.267);
    cmp_ok($result->get_statistic('entropy'), '==', 0.14);
    #is($result->get_statistic('dbletters'), 31984247);
    #is($result->get_statistic('dbentries'), 88780);
    #is($result->get_statistic('effective_hsplength'), 49);
    is($result->get_statistic('effectivespace'), '6.44279e+07');
    is($result->get_parameter('matrix'), 'BLOSUM62');
    is($result->get_parameter('gapopen'), 11);
    is($result->get_parameter('gapext'), 1);
    
    my $iter_count = 0;
    my @valid_hit_data = ( [ 'gi|16411294|emb|CAC99918.1|', 183, 'CAC99918', 16411294, '4.5377e-56', 209.92],
                   [ 'gi|16409584|emb|CAD00746.1|', 648, 'CAD00746', 16409584, '0.000286309', 37.7354],
                   [ 'gi|16411285|emb|CAC99909.1|', 209, 'CAC99909', 16411285, '0.107059', 29.261]);
    my @valid_iter_data = ( [ 16, 16, 0, 2, 14, 0, 0, 0, 0],
                [ 16, 8, 8, 0, 8, 0, 2, 0, 6]);
    
    while (my $iter = $result->next_iteration) {
        $iter_count++;
        my $di = shift @valid_iter_data;
        is($iter->number, $iter_count);
        is($iter->num_hits, shift @$di);
        is($iter->num_hits_new, shift @$di);
        is($iter->num_hits_old, shift @$di);
        is(scalar($iter->newhits_below_threshold), shift @$di);
        is(scalar($iter->newhits_not_below_threshold), shift @$di);
        is(scalar($iter->newhits_unclassified), shift @$di);
        is(scalar($iter->oldhits_below_threshold), shift @$di);
        is(scalar($iter->oldhits_newly_below_threshold), shift @$di);
        is(scalar($iter->oldhits_not_below_threshold), shift @$di);
        my $hit_count = 0;
        if ($iter_count == 1) {
            while( my $hit = $result->next_hit ) {
                my $d = shift @valid_hit_data;
                is($hit->name, shift @$d);
                is($hit->length, shift @$d);
                is($hit->accession, shift @$d);
                is($hit->ncbi_gi, shift @$d);
                is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
                is($hit->bits, shift @$d );
                if( $hit_count == 1 ) {
                    my $hsps_left = 1;
                    while( my $hsp = $hit->next_hsp ){
                        is($hsp->query->start, 4);
                        is($hsp->query->end, 155);
                        is($hsp->hit->start, 475);
                        is($hsp->hit->end, 617);
                        is($hsp->length('hsp'), 153);
                        is($hsp->start('hit'), $hsp->hit->start);
                        is($hsp->end('query'), $hsp->query->end);
                        is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
                        cmp_ok($hsp->evalue, '==', 0.000286309);
                        is($hsp->score, 86);
                        is($hsp->bits, 37.7354);
                        is(sprintf("%.1f",$hsp->percent_identity), 20.9);
                        is(sprintf("%.4f",$hsp->frac_identical('query')), 0.2105);
                        is(sprintf("%.3f",$hsp->frac_identical('hit')), 0.224);
                        is($hsp->gaps, 11);
                        $hsps_left--;
                    }
                    is($hsps_left, 0);
                }
                last if( $hit_count++ > @valid_hit_data );
            }
        }
    }
    is(@valid_hit_data, 0);
    is(@valid_iter_data, 0);
    is($iter_count, 2);
    
    $result = $searchio->next_result;
    is($result->database_name, 'AL591824.faa');
    is($result->database_entries, 2846);
    is($result->database_letters, 870878);
    is($result->algorithm, 'BLASTP');
    like($result->algorithm_version, qr/2\.2\.16/);
    is($result->query_name, 'gi|154350371|gb|ABS72450.1|');
    is($result->query_accession, 'ABS72450.1');
    is($result->query_gi, '154350371');
    is($result->query_length, 378);
    is($result->get_statistic('kappa') , 0.0450367);
    cmp_ok($result->get_statistic('lambda'), '==', 0.267);
    cmp_ok($result->get_statistic('entropy'), '==', 0.14);
    is($result->get_statistic('effectivespace'), '1.88702e+08');
    is($result->get_parameter('matrix'), 'BLOSUM62');
    is($result->get_parameter('gapopen'), 11);
    is($result->get_parameter('gapext'), 1);
    
    $iter_count = 0;
    
    @valid_hit_data = ( [ 'gi|16409361|emb|CAC98217.1|', 381, 'CAC98217', 16409361, '5.57178e-119', 420.239],
                   [ 'gi|16409959|emb|CAC98662.1|', 776, 'CAC98662', 16409959, '0.0242028', 32.7278],
                   [ 'gi|16410942|emb|CAC99591.1|', 382, 'CAC99591', 16410942, '0.340848', 28.8758]);
    @valid_iter_data = ( [ 11, 11, 0, 1, 10, 0, 0, 0, 0],
                [ 19, 11, 8, 0, 11, 0, 1, 0, 7]);
    
    while (my $iter = $result->next_iteration) {
        $iter_count++;
        my $di = shift @valid_iter_data;
        is($iter->number, $iter_count);
        is($iter->num_hits, shift @$di);
        is($iter->num_hits_new, shift @$di);
        is($iter->num_hits_old, shift @$di);
        is(scalar($iter->newhits_below_threshold), shift @$di);
        is(scalar($iter->newhits_not_below_threshold), shift @$di);
        is(scalar($iter->newhits_unclassified), shift @$di);
        is(scalar($iter->oldhits_below_threshold), shift @$di);
        is(scalar($iter->oldhits_newly_below_threshold), shift @$di);
        is(scalar($iter->oldhits_not_below_threshold), shift @$di);
        my $hit_count = 0;
        if ($iter_count == 1) {
            while( my $hit = $result->next_hit ) {
                my $d = shift @valid_hit_data;
                is($hit->name, shift @$d);
                is($hit->length, shift @$d);
                is($hit->accession, shift @$d);
                is($hit->ncbi_gi, shift @$d);
                is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
                is($hit->bits, shift @$d );
                if( $hit_count == 1 ) {
                    my $hsps_left = 1;
                    while( my $hsp = $hit->next_hsp ){
                        is($hsp->query->start, 63);
                        is($hsp->query->end, 181);
                        is($hsp->hit->start, 304);
                        is($hsp->hit->end, 432);
                        is($hsp->length('hsp'), 129);
                        is($hsp->start('hit'), $hsp->hit->start);
                        is($hsp->end('query'), $hsp->query->end);
                        is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
                        cmp_ok($hsp->evalue, '==', 0.0242028);
                        is($hsp->score, 73);
                        is($hsp->bits, 32.7278);
                        is(sprintf("%.1f",$hsp->percent_identity), '24.0');
                        is(sprintf("%.4f",$hsp->frac_identical('query')), '0.2605');
                        is(sprintf("%.3f",$hsp->frac_identical('hit')), '0.240');
                        is($hsp->gaps, 10);
                        $hsps_left--;
                    }
                    is($hsps_left, 0);
                }
                last if( $hit_count++ > @valid_hit_data );
            }
        }
    }
    is(@valid_hit_data, 0);
    is(@valid_iter_data, 0);
    is($iter_count, 2);
}

$searchio = Bio::SearchIO->new('-format' => 'blast',
				  '-file'   => test_input_file('ecolitst.bls'));

$result = $searchio->next_result;

is($result->database_name, 'ecoli.aa', 'database_name()');
is($result->database_entries, 4289);
is($result->database_letters, 1358990);

is($result->algorithm, 'BLASTP');
like($result->algorithm_version, qr/^2\.1\.3/);
like($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
is($result->query_accession, 'AAC73113.1');
is($result->query_gi, 1786183);
is($result->query_length, 820);
is($result->get_statistic('kappa'), '0.135');
is($result->get_statistic('kappa_gapped'), '0.0410');
is($result->get_statistic('lambda'), '0.319');
is($result->get_statistic('lambda_gapped'), '0.267');
is($result->get_statistic('entropy'), '0.383');
is($result->get_statistic('entropy_gapped'), '0.140');

is($result->get_statistic('dbletters'), 1358990);
is($result->get_statistic('dbentries'), 4289);
is($result->get_statistic('effective_hsplength'), 47);
is($result->get_statistic('effectivespace'), 894675611);
is($result->get_parameter('matrix'), 'BLOSUM62');
is($result->get_parameter('gapopen'), 11);
is($result->get_parameter('gapext'), 1);
is($result->get_statistic('S2'), '92');
is($result->get_statistic('S2_bits'), '40.0');
is($result->get_parameter('expect'), '1.0e-03');
is($result->get_statistic('num_extensions'), '82424');


my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 1567, 4058],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922', '1e-91', 332, 850],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994', '3e-47', 184, 467]);
my $count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->bits, shift @$d );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 820);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 820);
            is($hsp->length('hsp'), 820);
            is($hsp->start('hit'), $hsp->hit->start);
            is($hsp->end('query'), $hsp->query->end);
            is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            is($hsp->evalue, '0.0');
            is($hsp->score, 4058);
            is($hsp->bits,1567);	    	    
            is(sprintf("%.2f",$hsp->percent_identity), 98.29);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.9829);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9829);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('ecolitst.wublastp'));

$result = $searchio->next_result;

is($result->database_name, 'ecoli.aa');
is($result->database_letters, 1358990);
is($result->database_entries, 4289);
is($result->algorithm, 'BLASTP');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
like($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
is($result->query_accession, 'AAC73113.1');

is($result->query_length, 820);
is($result->query_gi, 1786183);
is($result->get_statistic('kappa'), 0.136);
is($result->get_statistic('lambda'), 0.319);
is($result->get_statistic('entropy'), 0.384);
is($result->get_statistic('dbletters'), 1358990);
is($result->get_statistic('dbentries'), 4289);
is($result->get_parameter('matrix'), 'BLOSUM62');
is($result->get_statistic('Frame+0_lambda_used'), '0.319');
is($result->get_statistic('Frame+0_kappa_used'), '0.136');
is($result->get_statistic('Frame+0_entropy_used'), '0.384');

is($result->get_statistic('Frame+0_lambda_computed'), '0.319');
is($result->get_statistic('Frame+0_kappa_computed'), '0.136');
is($result->get_statistic('Frame+0_entropy_computed'), '0.384');

is($result->get_statistic('Frame+0_lambda_gapped'), '0.244');
is($result->get_statistic('Frame+0_kappa_gapped'), '0.0300');
is($result->get_statistic('Frame+0_entropy_gapped'), '0.180');

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 4141],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922', '3.1e-86', 844],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994', '2.8e-47', 483]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    if ($count==1) {
        # Test HSP contig data returned by SearchUtils::tile_hsps()
        # Second hit has two hsps that overlap.
        my($qcontigs, $scontigs) = Bio::Search::SearchUtils::tile_hsps($hit);
        # Query contigs
        is($qcontigs->[0]->{'start'}, 5);
        is($qcontigs->[0]->{'stop'}, 812);
        is($qcontigs->[0]->{'iden'}, 250);
        is($qcontigs->[0]->{'cons'}, 413);
        # Subject contigs
        is($scontigs->[0]->{'start'}, 16);
        is($scontigs->[0]->{'stop'}, 805);
        is($scontigs->[0]->{'iden'}, 248);
        is($scontigs->[0]->{'cons'}, 410);
    }

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 820);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 820);
            is($hsp->length('hsp'), 820);
            
            is($hsp->evalue, '0.0');
            is($hsp->pvalue, '0.0');
            is($hsp->score, 4141);
            is($hsp->bits,1462.8);	    	    
            is($hsp->percent_identity, 100);
            is($hsp->frac_identical('query'), 1.00);
            is($hsp->frac_identical('hit'), 1.00);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

# test that add hit really works properly for BLAST objects
# bug 1611
my @hits = $result->hits;
$result->add_hit($hits[0]);
is($result->num_hits, @hits + 1);

# test WU-BLAST -noseqs option
$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('ecolitst.noseqs.wublastp'));

$result = $searchio->next_result;

is($result->database_name, 'ecoli.aa');
is($result->database_letters, 1358990);
is($result->database_entries, 4289);
is($result->algorithm, 'BLASTP');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
like($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
is($result->query_accession, 'AAC73113.1');
is($result->query_gi, 1786183);

is($result->query_length, 820);
is($result->get_statistic('kappa'), 0.135);
is($result->get_statistic('lambda'), 0.319);
is($result->get_statistic('entropy'), 0.384);
is($result->get_statistic('dbletters'), 1358990);
is($result->get_statistic('dbentries'), 4289);
is($result->get_parameter('matrix'), 'BLOSUM62');
is($result->get_statistic('Frame+0_lambda_used'), '0.319');
is($result->get_statistic('Frame+0_kappa_used'), '0.135');
is($result->get_statistic('Frame+0_entropy_used'), '0.384');

is($result->get_statistic('Frame+0_lambda_computed'), '0.319');
is($result->get_statistic('Frame+0_kappa_computed'), '0.135');
is($result->get_statistic('Frame+0_entropy_computed'), '0.384');

is($result->get_statistic('Frame+0_lambda_gapped'), '0.244');
is($result->get_statistic('Frame+0_kappa_gapped'), '0.0300');
is($result->get_statistic('Frame+0_entropy_gapped'), '0.180');

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 4141],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922', '6.6e-93', 907],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994', '2.8e-47', 483]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 820);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 820);
            is($hsp->length('hsp'), 820);
            
            is($hsp->evalue , '0.');
            is($hsp->pvalue , '0.');
            is($hsp->score, 4141);
            is($hsp->bits,1462.8);	    	    
            is($hsp->percent_identity, 100);
            is($hsp->frac_identical('query'), 1.00);
            is($hsp->frac_identical('hit'), 1.00);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

# test tblastx 
$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('HUMBETGLOA.tblastx'));

$result = $searchio->next_result;
is($result->database_name, 'ecoli.nt');
is($result->database_letters, 4662239);
is($result->database_entries, 400);
is($result->algorithm, 'TBLASTX');
like($result->algorithm_version, qr/^2\.1\.2/);
is($result->query_name, 'HUMBETGLOA');
is($result->query_description, 'Human haplotype C4 beta-globin gene, complete cds.');
is($result->query_length, 3002);
is($result->get_statistic('kappa'), 0.135);
is($result->get_statistic('lambda'), 0.318);
is($result->get_statistic('entropy'), 0.401);
is($result->get_statistic('dbletters'), 4662239);
is($result->get_statistic('dbentries'), 400);
is($result->get_statistic('T'), 13);
is($result->get_statistic('X1'), 16);
is($result->get_statistic('X1_bits'), 7.3);
is($result->get_statistic('X2'), 0);
is($result->get_statistic('X2_bits'), '0.0');
is($result->get_statistic('S1'), 41);
is($result->get_statistic('S1_bits'), 21.7);
is($result->get_statistic('S2'), 53);
is($result->get_statistic('S2_bits'), 27.2);

is($result->get_statistic('decayconst'), 0.1);

is($result->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AE000479.1|AE000479', 10934, 'AE000479', '0.13', 33.6, 67],
	   [ 'gb|AE000302.1|AE000302', 10264, 'AE000302', '0.61', 31.3, 62],
	   [ 'gb|AE000277.1|AE000277', 11653, 'AE000277', '0.84', 30.8, 61]);
$count = 0;

while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->significance, shift @$d );
    is($hit->bits, shift @$d );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1057);
            is($hsp->query->end, 1134);
            is($hsp->query->strand, 1);
            is($hsp->strand('query'), $hsp->query->strand);
            is($hsp->hit->end, 5893);
            is($hsp->hit->start, 5816);
            is($hsp->hit->strand, -1);
            is($hsp->strand('sbjct'), $hsp->subject->strand);
            is($hsp->length('hsp'), 26);
            
            is($hsp->evalue , 0.13);
            is($hsp->score, 67);
            is($hsp->bits,33.6);
            is(sprintf("%.2f",$hsp->percent_identity), 42.31);
            is(sprintf("%.4f",$hsp->frac_identical('query')), '0.4231');
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.4231');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 1);
            is($hsp->gaps, 0);	    
            is($hsp->query_string, 'SAYWSIFPPLGCWWSTLGPRGSLSPL');
            is($hsp->hit_string, 'AAVWALFPPVGSQWGCLASQWRTSPL');
            is($hsp->homology_string, '+A W++FPP+G  W  L  +   SPL');
            # changed to reflect positional ambiguities, note extra flag
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '1063-1065 1090-1095 1099-1104 1108-1113 1117-1125');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '5825-5833 5837-5842 5846-5851 5855-5860 5885-5887');            
            is($hsp->ambiguous_seq_inds, 'query/subject');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('HUMBETGLOA.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/dros_clones.2.5/);
is($result->database_letters, 112936249);
is($result->database_entries, 657);
is($result->algorithm, 'FASTN');
is($result->algorithm_version, '3.3t08');
is($result->query_name, "HUMBETGLOA");
is($result->query_description, "Human haplotype C4 beta-globin gene, complete cds.");
is($result->query_length, 3002);
is($result->get_parameter('gapopen'), -16);
is($result->get_parameter('gapext'), -4);
is($result->get_parameter('ktup'), 6);

is($result->get_statistic('lambda'), 0.0823);
is($result->get_statistic('dbletters'), 112936249);
is($result->get_statistic('dbentries'), 657);

@valid = ( [ 'BACR21I23', 73982, 'BACR21I23', '0.017', 44.2],
	   [ 'BACR40P19', 73982, 'BACR40P19', '0.017', 44.2],
	   [ 'BACR30L17', 32481, 'BACR30L17', '0.018', 44.1]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->significance, shift @$d );
    is($hit->raw_score, shift @$d );
    is($hit->rank, $count + 1);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 31);
            is($hsp->query->end, 289);
            is($hsp->query->strand, -1);
            is($hsp->hit->end, 65167);
            is($hsp->hit->start, 64902);
            is($hsp->hit->strand, 1);
            is($hsp->length('hsp'), 267);	    
            is($hsp->evalue , 0.017);
            is($hsp->score, 134.5);
            is($hsp->bits,44.2);
            is(sprintf("%.2f",$hsp->percent_identity), '57.30');
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5907); 
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5752);
	    # these are really UNGAPPED values not CONSERVED
	    # otherwise ident and conserved would be identical for
	    # nucleotide alignments
	    is(sprintf("%.4f",$hsp->frac_conserved('total')), 0.5955); 
	    is(sprintf("%.4f",$hsp->frac_conserved('query')), 0.6139); 
	    is(sprintf("%.4f",$hsp->frac_conserved('hit')), 0.5977); 
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps, 159);
            is($hsp->gaps('query'), 8);
            is($hsp->gaps('hit'),1);
            is($hsp->query_string, 'GATTAAAACCTTCTGGTAAGAAAAGAAAAAATATATATATATATATATGTGTATATGTACACACATACATATACATATATATGCATTCATTTGTTGTTGTTTTTCTTAATTTGCTCATGCATGCTA----ATAAATTATGTCTAAAAATAGAAT---AAATACAAATCAATGTGCTCTGTGCATTA-GTTACTTATTAGGTTTTGGGAAACAAGAGGTAAAAAACTAGAGACCTCTTAATGCAGTCAAAAATACAAATAAATAAAAAGTCACTTACAACCCAAAGTGTGACTATCAATGGGGTAATCAGTGGTGTCAAATAGGAGGT');
            is($hsp->hit_string, 'GATGTCCTTGGTGGATTATGGTGTTAGGGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATATAAAATATAATATAAAATATAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATAT-AATATAAAATATAAAATAAAATATAATATAAAATATAATATAAAATATAATATAAAATATAATATAAAATA');
            is($hsp->homology_string, '                              :::::::::::::::::: : ::::: :: : : ::: ::::: ::::::::  ::  :: : :   : : : : :  ::    : :: ::   ::    : ::: :::     :::::: :::   ::::: ::  :::  :    :    : ::   :::  : ::   : :   : : :: :   :: : : :: : :       ::  : : ::: ::: ::  ::::: ::: : :  :: ::   ::: : : : ::: ::   '.' 'x60);
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '33 37 39 41 43 47-49 52 55 56 58 60 64 70 71 74 78 82 84 86 87 90-96 98 100 103 105 107 110-112 114 117 119 121-123 125 127-129 132 134 135 139-141 143 145-148 150-153 155 156 160 161 164 170 173 180-184 188 192 194 196-198 201 204 206-209 212 213 215 217 219 221 223-225 227 229 232 233 236 237 246 252 256 258 260 263 269 271');
            is(join(' ', $hsp->seq_inds('query', 'conserved',1)), '31 32 34-36 38 40 42 44-46 50 51 53 54 57 59 61-63 65-69 72 73 75-77 79-81 83 85 88 89 97 99 101 102 104 106 108 109 113 115 116 118 120 124 126 130 131 133 136-138 141 142 144 149 154 157-159 162 163 165-172 174-179 185-187 189-191 193-195 199 200 202 203 205 210 211 214 216 218 220 222 226 228 230 231 234 235 238-245 247-251 253-255 257 259 261 262 264-268 270 272-289');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '64920 64922 64928 64931 64933 64935 64939 64945 64954 64955 64958 64959 64962 64964 64966-64968 64970 64972 64974 64976 64978 64979 64982-64985 64987 64990 64993-64995 64998-65001 65003 65007 65011-65015 65022 65026-65028 65034 65037 65038 65042 65043 65045-65048 65050-65053 65055 65058-65060 65064 65065 65067 65070-65072 65074 65076-65078 65080 65082 65085 65087-65089 65092 65094 65096 65099 65101 65103-65109 65112 65113 65115 65117 65121 65125 65128 65129 65135 65139 65141 65143 65144 65147 65150-65152 65156 65158 65159 65161 65165');
            is(join(' ', $hsp->seq_inds('hit', 'conserved',1)), '64902-64919 64921 64923-64927 64929 64930 64932 64934 64936-64938 64940-64944 64946-64953 64956 64957 64960 64961 64963 64965 64969 64971 64973 64975 64977 64980 64981 64986 64988 64989 64991 64992 64996 64997 65002 65004-65006 65008-65010 65016-65021 65023-65025 65029-65033 65035 65036 65039-65041 65044 65049 65054 65056 65057 65061-65063 65066 65068 65069 65073 65075 65079 65081 65083 65084 65086 65090 65091 65093 65095 65097 65098 65100 65102 65110 65111 65114 65116 65118-65120 65122-65124 65126 65127 65130-65134 65136-65138 65140 65142 65145 65146 65148 65149 65153-65155 65157 65159 65160 65162-65164 65166 65167');            
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '141 170 194');
            is($hsp->ambiguous_seq_inds, '');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), '59.30');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('cysprot1.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/ecoli.aa/);
is($result->database_letters, 1358987);
is($result->database_entries, 4289);
is($result->algorithm, 'FASTP');
is($result->algorithm_version, '3.3t08');
is($result->query_name, 'CYS1_DICDI');
is($result->query_length, 343);
is($result->get_parameter('gapopen'), -12);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);

is($result->get_statistic('lambda'), 0.1456);
is($result->get_statistic('dbletters'), 1358987);
is($result->get_statistic('dbentries'), 4289);


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309', 1787478, 1.2, 29.2],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148', 1790635, 2.1, 27.4],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494', 1786590, 2.1, 25.9]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->ncbi_gi, shift @$d);
    is($hit->significance, shift @$d );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 125);
            is($hsp->query->end, 305);
            is($hsp->query->strand, 0);
            is($hsp->hit->start, 2);
            is($hsp->hit->end, 181);
            is($hsp->hit->strand, 0);
            is($hsp->length('hsp'), 188);	    
            is($hsp->evalue , 1.2);
            is($hsp->score, 109.2);
            is($hsp->bits,29.2);
            is(sprintf("%.2f",$hsp->percent_identity), 23.94);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.2486);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.2500');
            is(sprintf("%.4f",$hsp->frac_conserved('query')), '0.2707');
            is(sprintf("%.4f",$hsp->frac_conserved('hit')), '0.2722');
	    # there is slight rounding different here so file says 26.012%
	    # but with the rounding this ends up as 0.2606
            is(sprintf("%.4f",$hsp->frac_conserved('total')), '0.2606');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 7);
            is($hsp->gaps, 49);	    
            is($hsp->query_string, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
            is($hsp->hit_string, (' 'x29).'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
            is($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                      '.' 'x60);
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 26.01);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

is($result->hits, 8);
$searchio = Bio::SearchIO->new(-format => 'fasta',
				 -file   => test_input_file('cysprot_vs_gadfly.FASTA'));
$result = $searchio->next_result;
like($result->database_name, qr/gadflypep2/);
is($result->database_letters, 7177762);
is($result->database_entries, 14334);
is($result->algorithm, 'FASTP');
is($result->algorithm_version, '3.3t08');
is($result->query_name, 'cysprot.fa');
is($result->query_length, 2385);
is($result->get_parameter('gapopen'), -12);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');

is($result->get_statistic('lambda'), 0.1397);
is($result->get_statistic('dbletters'), 7177762 );
is($result->get_statistic('dbentries'), 14334);


@valid = ( [ 'Cp1|FBgn0013770|pp-CT20780|FBan0006692', 341, 
	     'FBan0006692', '3.1e-59', 227.8],
	   [ 'CG11459|FBgn0037396|pp-CT28891|FBan0011459', 336, 
	     'FBan0011459', '6.4e-41',  166.9],
	   [ 'CG4847|FBgn0034229|pp-CT15577|FBan0004847', 390, 
	     'FBan0004847',  '2.5e-40', 165.2]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1373);
            is($hsp->query->end, 1706);
            is($hsp->query->strand, 0);
            is($hsp->hit->start, 5);
            is($hsp->hit->end, 341);
            is($hsp->hit->strand, 0);
            is($hsp->length('hsp'), 345);	    
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'3.1e-59') );
            is($hsp->score, 1170.6);
            is($hsp->bits,227.8);
            is(sprintf("%.2f",$hsp->percent_identity), 53.04);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5479);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.5430');
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 11);
            is($hsp->gaps, 194);
            is($hsp->hit_string, (' 'x26).'MRTAVLLPLLAL----LAVAQA-VSFADVVMEEWHTFKLEHRKNYQDETEERFRLKIFNENKHKIAKHNQRFAEGKVSFKLAVNKYADLLHHEFRQLMNGFNYTLHKQLRAADESFKGVTFISPAHVTLPKSVDWRTKGAVTAVKDQGHCGSCWAFSSTGALEGQHFRKSGVLVSLSEQNLVDCSTKYGNNGCNGGLMDNAFRYIKDNGGIDTEKSYPYEAIDDSCHFNKGTVGATDRGFTDIPQGDEKKMAEAVATVGPVSVAIDASHESFQFYSEGVYNEPQCDAQNLDHGVLVVGFGTDESGED---YWLVKNSWGTTWGDKGFIKMLRNKENQCGIASASSYPLV');
            is($hsp->query_string, 'SNWGNNGYFLIERGKNMCGLAACASYPIPQVMNPTLILAAFCLGIASATLTFDHSLEAQWTKWKAMHNRLY-GMNEEGWRRAVWEKNMKMIELHNQEYREGKHSFTMAMNAFGDMTSEEFRQVMNGFQ---NRKPR------KGKVFQEPLFYEAPRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQDNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPK-QEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCSSEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTVMTPLLLLAVLCLGTALATPKFDQTFNAQWHQWKSTHRRLYGTNEE');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 56.13);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);
is($result->hits, 21);

# test on TFASTXY
$searchio = Bio::SearchIO->new(-format => 'fasta',
			      -file   => test_input_file('5X_1895.FASTXY'));
$result = $searchio->next_result;
like($result->database_name, qr/yeast_nrpep.fasta/);
is($result->database_letters, 4215311);
is($result->database_entries, 9190);
is($result->algorithm, 'FASTY');
is($result->algorithm_version, '3.4t07');
is($result->query_name, '5X_1895.fa');
is($result->query_length, 7972);
is($result->get_parameter('gapopen'), -14);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');

is($result->get_statistic('lambda'), 0.1711);
is($result->get_statistic('dbletters'), 4215311);
is($result->get_statistic('dbentries'), 9190);


@valid = ( [ 'NR_SC:SW-YNN2_YEAST', 1056, 'NR_SC:SW-YNN2_YEAST','1.6e-154', '547.0'],
	   [ 'NR_SC:SW-MPCP_YEAST', 311, 'NR_SC:SW-MPCP_YEAST', '1.3e-25', 117.1],
	   [ 'NR_SC:SW-YEO3_YEAST', 300, 'NR_SC:SW-YEO3_YEAST', '5.7e-05', 48.5]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 2180);
            is($hsp->query->end, 5623);
            is($hsp->query->strand, 1);
            is($hsp->hit->start, 3);
            is($hsp->hit->end, 1053);
            is($hsp->hit->strand, 0);
            is($hsp->length('hsp'), 1165);
            
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'1.6e-154'));
            is($hsp->score, 2877.6);
            is($hsp->bits,'547.0');
            is(sprintf("%.2f",$hsp->percent_identity), 51.67);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.5244);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5728);
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps, 678);	    
            is($hsp->query_string, 'RKQLDPRIPALINNGVKANHRSFFVMVGDKGRDQVCPGMQAAMRFD*HRCR/LVNLHFLLSQARVSSRPSVLWCYKKD-LGFTT*VAASENLQQTIYFRPIATSHRKKREAKIKRDVKRGIRDANEQDPFELFVTVTDIRYTYYKDSAKILGQTFGMLVLQDYEAITPNLLARTIETVEGGGIVVLLLKTMSSLKQLYAMAM/DKL*CRDGVE*SDFS*LLI*DVHSRYRTDAHQFVQPRFNERFILSLGSNPDCLVLDDELNVLPLSKGKDIQIGKAGEEDDRGRKRKAEELKEMKENLEGVDIVGSLAKLAKTVDQAKAILTFVEAISEKNLSSTVALTAGRGRGKSAALGLAIGAALAHDYSNIFVTSPDPENLKTLFEFVFKALDALGYEEHIDYDVVQSTNPDFKKAIVRVNIFRGHRQTIQYISPEDSHVLGQAELVIIDEAAAIPLPLVRKLIGPYLVFMASTINGYEGTGRSLSIKLIQQLREQTRPSITKDSENAAASSAGSSSKAAAAGRSGAGLVRSLREIKLDEPIRYSPGDNVEKWLNNLLCLDATIVSK---SIQGCPHPSKCELYYVNRDTLFSYHPASEVFLQRMMALYVASHYKNSPNDLQMLSDAPAHHLFVLLPPIDEND-NTLPDPLVVLQVALEGNISREAILKEMAQSGMRSSGDMIPWIISTQFQDNDFATLSGARVVRIATHPDYARMGYGSRAMEALESFYNGTSYNFDDVPVDMGESFAD\VPRSDL*VTSFIPFPQNRTSTECVSQNANLQNDTIAIRDPSRMPPLLQRLSERKPETLDYLGVSFGLTRDLLRFWKKGGFTPLYASQKENALTGEYTFVMLKVLASAGGGGEWLGAFAQGMSCLLLQDEVHMGND*RL*TDFRQRFMNLLSYEAFKKFDASIALSILESTVPRNSPSPAP----KLLTNTELSSLLTPFDIKRLESYADSMLDYHVVLDLVPTIASLFFGKRLETS--LPPAQQAILLALGLQRKNVEALENELGITSTQTLALFGKVLRKMTKSLEDIRKASIASELP-----AEPTLAGRSANGSNKFVALQQTIEQDLADSAVQLNGEDDDASKKEQRELLNTLNMEEFAI-DQGGDWTEAEKQVERLASGKGGTRLSSTVSVKVDKLDD\AKRRRRRARMRVPRMRRR');
            is($hsp->hit_string, 'KKAIDSRIPSLIRNGVQTKQRSIFVIVGDRARNQ------------------LPNLHYLMMSADLKMNKSVLWAYKKKLLGFT--------------------SHRKKRENKIKKEIKRGTREVNEMDPFESFISNQNIRYVYYKESEKILGNTYGMCILQDFEALTPNLLARTIETVEGGGIVVILLKSMSSLKQLYTMTM-D--------------------VHARYRTEAHGDVVARFNERFILSLGSNPNCLVVDDELNVLPLSGAKNVKPLPPKEDDELPPKQL--ELQELKESLEDVQPAGSLVSLSKTVNQAHAILSFIDAISEKTLNFTVALTAGRGRGKSAALGISIAAAVSHGYSNIFVTSPSPENLKTLFEFIFKGFDALGYQEHIDYDIIQSTNPDFNKAIVRVDIKRDHRQTIQYIVPQDHQVLGQAELVVIDEAAAIPLPIVKNLLGPYLVFMASTINGYEGTGRSLSLKLIQQLRNQNNTSGRESTQTAVVSRDNKEKDSHLHSQS-----RQLREISLDEPIRYAPGDPIEKWLNKLLCLDVTLIKNPRFATRGTPHPSQCNLFVVNRDTLFSYHPVSENFLEKMMALYVSSHYKNSPNDLQLMSDAPAHKLFVLLPPIDPKDGGRIPDPLCVIQIALEGEISKESVRNSLSR-GQRAGGDLIPWLISQQFQDEEFASLSGARIVRIATNPEYASMGYGSRAIELLRDYFEGKF-------TDMSE---D-VRPKDYSI--------KRVSDKELAKT-NLLKDDVKLRDAKTLPPLLLKLSEQPPHYLHYLGVSYGLTQSLHKFWKNNSFVPVYLRQTANDLTGEHTCVMLNVLE--GRESNWLVEFAK---------------------DFRKRFLSLLSYD-FHKFTAVQALSVIESSKKAQDLSDDEKHDNKELTRTHLDDIFSPFDLKRLDSYSNNLLDYHVIGDMIPMLALLYFGDKMGDSVKLSSVQSAILLAIGLQRKNIDTIAKELNLPSNQTIAMFAKIMRKMSQYFRQLLSQSIEETLPNIKDDAIAEMDGEEIKNYNAAEALDQ-MEEDLEEAG----SEAVQAMREKQKELINSLNLDKYAINDNSEEWAESQKSLEIAAKAKGVVSLKTGKKRTTEKAED-IYRQEMKA-MKKPRKSKK');
            is($hsp->homology_string, '.: .: :::.:: :::....::.::.:::..:.:                  : :::.:. .: ..   :::: :::  ::::                    ::::::: :::...::: :..::.:::: :..  .:::.:::.: ::::.:.:: .:::.::.:::::::::::::::::::.:::.::::::::.:.: :                    ::.::::.::  :  ::::::::::::::.:::.:::::::::: .:...     :.:.   :.   ::.:.::.:: :. .:::..:.:::.::.:::.:..:::::.:. :::::::::::::::::..:.::..: :::::::::.::::::::::.::..:::::.::::::..:::::::.::::::.: : :::::::: :.: .::::::::.::::::::::.:..:.::::::::::::::::::::::.:::::::.:.  :  .....:..:  .. . .   ..:     :.::::.:::::::.::: .:::::.:::::.:....   . .: ::::.:.:. :::::::::::.:: ::..::::::.:::::::::::..::::::.::::::::: .: . .:::: :.:.::::.::.:.. . ... :.:..::.:::.:: ::::..::.:::::.:::::.:.:: :::::::.: :.....:         .::.:   : :  .:  .        .:.: . .... :: .: . .:: . .:::: .:::. :. : :::::.:::..: .:::...:.:.:  :  : ::::.: :::.::   :  ..::  ::.                     :::.::..::::. :.:: :  :::..::.   .. :       : :: :.:.....:::.:::.::....:::::. :..: .: :.:: ..  :  :  .:.:::::.::::::.... .::.. :.::.:.:.:..:::.. .... . ::   ::     :   . :.  .. :   ::.: .:.:: ...    .:  .: ...:.::.:.::....:: :.. .:.:..:..:  :..:: . :..  .  ..: .:   :.. .: :. ::  ..');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            is(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity()),
               '51.77');
            is(sprintf("%.2f",$hsp->get_aln->average_percentage_identity()),
               '58.41');
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '2186-2188 2195-2197 2216-2218 2282-2335 2339-2341 2360-2362 2369-2371 2378-2386 2399-2401 2411-2416 2426-2485 2507-2509 2537-2539 2570-2572 2582-2587 2618-2620 2648-2650 2783-2785 2789-2848 2879-2884 2888-2893 2981-2983 2999-3013 3026-3034 3041-3049 3080-3082 3089-3091 3182-3184 3263-3265 3431-3433 3437-3439 3464-3466 3476-3478 3656-3661 3665-3670 3698-3703 3710-3712 3716-3718 3722-3730 3740-3754 3809-3811 3866-3871 3878-3880 3908-3910 3953-3955 4076-4078 4085-4090 4106-4108 4154-4156 4160-4162 4172-4174 4217-4219 4295-4297 4325-4327 4349-4375 4391-4399 4403-4405 4409-4414 4421-4426 4430-4453 4466-4468 4472-4474 4487-4489 4496-4498 4505-4507 4511-4513 4523-4525 4529-4531 4547-4549 4565-4567 4574-4576 4580-4582 4619-4621 4658-4663 4667-4672 4676-4678 4697-4699 4718-4726 4730-4735 4748-4753 4763-4825 4865-4867 4880-4882 4886-4891 4916-4924 4931-4933 4937-4951 4958-4960 5045-5047 5060-5062 5069-5071 5084-5086 5093-5098 5102-5110 5168-5170 5186-5188 5240-5242 5255-5257 5261-5263 5270-5278 5285-5296 5300-5302 5309-5314 5321-5323 5327-5335 5348-5350 5366-5368 5378-5389 5396-5401 5408-5410 5465-5467 5474-5476 5507-5512 5528-5530 5534-5536 5546-5551 5555-5560 5570-5572 5579-5587 5597-5599 5606-5608 5615-5617');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '9 12 19 40 42 49 52 55-57 62 66 67 71 79 89 100 104 105 116 126 170 171 182 183 185 186 216 222-226 231-233 236 247 250 281 308 364 366 375 379 439 440 442 443 453 454 457 459 461-463 466 485 504-506 508 511 521 536 577 580 582 588 604 606 609 624 650 660 668 669 674 675 677 678 681-683 688 690 694 697 700 702 706 708 714 720 723 725 738 751 752 754 755 757 764 771 773 774 779 780 783 796 801 803 804 813-815 818 820-826 828 831 860 865 868 873 876 877 879 880 882 883 903 909 927 932 934 937-939 942-946 948-950 952 955 956 959 961-963 967 973 976 979 980 983 1002 1006 1017 1018 1024 1026 1030 1031 1033 1034 1038 1040-1042 1046 1048 1051');
            is(join(' ', $hsp->seq_inds('query', 'gap',1)), '2422 3871 4090 4948 5104 5287 5467');
            is(join(' ', $hsp->seq_inds('hit', 'gap',1)), '40 71 170 171 236 466 609 669 674 675 683 694 771 783 796 967 976 1040 1048');
            is($hsp->ambiguous_seq_inds, 'query');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);
is($result->hits, 58);

# test FASTA v35.04, params encoding changed 
# test on TFASTXY
$searchio = Bio::SearchIO->new(-format => 'fasta',
			      -file   => test_input_file('BOSS_DROME.FASTP_v35_04'));
$result = $searchio->next_result;
like($result->database_name, qr/wormpep190/);
is($result->database_letters, 10449259);
is($result->database_entries, 23771);
is($result->algorithm, 'FASTA');
is($result->algorithm_version, '35.04');
is($result->query_name, 'BOSS_DROME');
is($result->query_length, 896);
is($result->get_parameter('gapopen'), -10);
is($result->get_parameter('gapext'), -2);
is($result->get_parameter('ktup'), 2);
is($result->get_parameter('matrix'), 'BL50');
is($result->get_parameter('wordsize'), 16);
is($result->get_parameter('filter'), '15:-5');

is($result->get_statistic('lambda'), 0.122629);
is($result->get_statistic('dbletters'), 10449259);
is($result->get_statistic('dbentries'), 23771);
is($result->get_statistic('effectivespace'),23771);

# test for MarkW bug in blastN

$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('a_thaliana.blastn'));


$result = $searchio->next_result;
is($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS, GSS,or phase 0, 1 or 2 HTGS sequences) ');
is($result->database_letters, 4677375331);
is($result->database_entries, 1083200);
is($result->algorithm, 'BLASTN');
like($result->algorithm_version, qr/^2\.2\.1/);
is($result->query_name, '');
is($result->query_length, 60);
is($result->get_parameter('gapopen'), 5);
is($result->get_parameter('gapext'), 2);
is($result->get_parameter('ktup'), undef);

is($result->get_statistic('lambda'), 1.37);
is($result->get_statistic('kappa'), 0.711);
is($result->get_statistic('entropy'),1.31 );
is($result->get_statistic('T'), 0);
is($result->get_statistic('A'), 30);
is($result->get_statistic('X1'), '6');
is($result->get_statistic('X1_bits'), 11.9);
is($result->get_statistic('X2'), 15);
is($result->get_statistic('X2_bits'), 29.7);
is($result->get_statistic('S1'), 12);
is($result->get_statistic('S1_bits'), 24.3);
is($result->get_statistic('S2'), 17);
is($result->get_statistic('S2_bits'), 34.2);

is($result->get_statistic('dbentries'), 1083200);

@valid = ( [ 'gb|AY052359.1|', 2826, 'AY052359', '3e-18', 95.6, 48, 1, 60, 
	     '1.0000'],
	   [ 'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 95.6, 48, 1, 60, 
	     '1.0000' ],
	   [ 'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.04', 42.1, 21, 35, 55, 
	     '0.3500']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->bits, shift @$d );
    is($hit->raw_score, shift @$d );
    is($hit->start, shift @$d);
    is($hit->end,shift @$d);    
    is(sprintf("%.4f",$hit->frac_aligned_query), shift @$d);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 60);
            is($hsp->query->strand, 1);
            is($hsp->hit->start, 154);
            is($hsp->hit->end, 212);
            is($hsp->hit->strand, 1);
            is($hsp->length('hsp'), 60);	    
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'3e-18'));
            is($hsp->score, 48);
            is($hsp->bits,95.6);
            is(sprintf("%.2f",$hsp->percent_identity), 96.67);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.9667);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9831);
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 0);
            is($hsp->gaps('hit'), 1);
            is($hsp->gaps, 1);	    
            is($hsp->query_string, 'aggaatgctgtttaattggaatcgtacaatggagaatttgacggaaatagaatcaacgat');
            is($hsp->hit_string, 'aggaatgctgtttaattggaatca-acaatggagaatttgacggaaatagaatcaacgat');
            is($hsp->homology_string, '|||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||');
			my $aln = $hsp->get_aln;
            is(sprintf("%.2f", $aln->overall_percentage_identity), 96.67);
            is(sprintf("%.2f",$aln->percentage_identity), 98.31);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
} 
is(@valid, 0);

#WU-BlastX test

$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('dnaEbsub_ecoli.wublastx'));

$result = $searchio->next_result;
is($result->database_name, 'ecoli.aa');
is($result->database_letters, 1358990);
is($result->database_entries, 4289);
is($result->algorithm, 'BLASTX');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
is($result->query_name, 'gi|142864|gb|M10040.1|BACDNAE');
is($result->query_description, 'B.subtilis dnaE gene encoding DNA primase, complete cds');
is($result->query_accession, 'M10040.1');
is($result->query_gi, 142864);
is($result->query_length, 2001);
is($result->get_parameter('matrix'), 'blosum62');

is($result->get_statistic('lambda'), 0.318);
is($result->get_statistic('kappa'), 0.135);
is($result->get_statistic('entropy'),0.401 );

is($result->get_statistic('dbentries'), 4289);

@valid = ( [ 'gi|1789447|gb|AAC76102.1|', 581, 'AAC76102', '1.1e-74', 671]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );
    is(sprintf("%.4f",$hit->frac_identical('query')), '0.3640');
    is(sprintf("%.4f",$hit->frac_identical('hit')), '0.3660');
    is(sprintf("%.4f",$hit->frac_conserved('query')), '0.5370');
    is(sprintf("%.4f",$hit->frac_conserved('hit')), '0.5400');
    is(sprintf("%.4f",$hit->frac_aligned_query), '0.6200');
    is(sprintf("%.4f",$hit->frac_aligned_hit), '0.7100');

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 21);
            is($hsp->query->end, 1265);
            is($hsp->query->strand, 1);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 413);
            is($hsp->hit->strand, 0);
            is($hsp->length('hsp'), 421);	    
            is(sprintf("%g",$hsp->evalue), sprintf("%g",'1.1e-74'));
            is(sprintf("%g",$hsp->pvalue), sprintf("%g",'1.1e-74'));
            is($hsp->score,671);
            is($hsp->bits,265.8);
            is(sprintf("%.2f",$hsp->percent_identity), 35.87);
            
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.3639);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.3656);
            is(sprintf("%.4f",$hsp->frac_conserved('query')), 0.5373);
            is(sprintf("%.2f",$hsp->frac_conserved('hit')), 0.54);
            
            is(sprintf("%.4f",$hsp->frac_identical('hsp')), 0.3587);
            is(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.5297);
            
            is($hsp->query->frame(), 2);
            is($hsp->hit->frame(), 0);
            is($hsp->gaps('query'), 6);
            is($hsp->gaps('hit'), 8);
            is($hsp->gaps, 14);	    
            is($hsp->query_string, 'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ');
            is($hsp->hit_string, 'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ');
            is($hsp->homology_string, 'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q');
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '24-29 39-47 54-56 60-71 90-98 129-137 150-152 156-158 180-182 192-194 219-221 228-236 243-251 255-263 267-269 279-284 291-296 300-302 309-311 315-317 321-332 342-347 351-362 366-368 372-374 378-383 387-389 393-398 405-413 417-440 444-449 456-461 468-470 474-476 486-491 495-497 510-518 525-527 531-533 537-557 561-569 573-578 594-599 603-605 609-614 618-620 633-635 654-656 660-665 669-671 678-680 684-686 693-695 705-710 738-740 753-755 759-761 768-773 786-797 801-806 810-812 819-821 831-833 840-860 864-869 894-896 900-902 921-923 927-938 945-947 957-959 972-974 981-986 993-995 999-1013 1017-1019 1029-1037 1050-1052 1062-1067 1077-1079 1083-1085 1089-1091 1098-1103 1107-1109 1113-1115 1122-1124 1128-1130 1137-1163 1173-1184 1188-1208 1212-1217 1224-1226 1230-1232 1236-1244 1248-1250');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '2 3 7-9 12 14-17 24-26 37-39 44 46 54 58 67 70-72 75-77 79-81 83 87 88 91 92 94 97 99 100 104 106-108 110-113 115 117 119 120 122 124 125 128-138 140 141 144 145 148 150 154 155 157 162-164 167 169 171-177 179-181 183 184 190 191 193 195-197 202 209 211 212 214 217 219 222 226 227 237 242 244 247 248 253-256 258 259 261 264 268 271-277 279 280 289 291 298 300-303 306 310 315 318 319 322 324-331 333 337-339 344 348 349 353 355 357 360-362 364 367 369 372-380 384-387 389-395 397 398 401 403 405-408');
            is($hsp->ambiguous_seq_inds, 'query');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

#Trickier WU-Blast
$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('tricky.wublast'));
$result = $searchio->next_result;
my $hits_left = 1;
while (my $hit = $result->next_hit) {
	# frac_aligned_hit used to be over 1, frac_identical & frac_conserved are still too wrong
	TODO: {
        local $TODO = 'frac_identical & frac_conserved are still too wrong';
        cmp_ok sprintf("%.3f",$hit->frac_identical), '>', 0.9;
        cmp_ok sprintf("%.3f",$hit->frac_conserved), '<=', 1;
    }
    is(sprintf("%.2f",$hit->frac_aligned_query), '0.92');
    is(sprintf("%.2f",$hit->frac_aligned_hit), '0.91');
    $hits_left--;
}
is($hits_left, 0);

# More frac_ method testing, this time on ncbi blastn
$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('frac_problems.blast'));
my @expected = ("1.000", "0.943");
while (my $result = $searchio->next_result) {
    my $hit = $result->next_hit;
    is($hit->frac_identical, shift @expected);
}
is(@expected, 0);

# And even more: frac_aligned_query should never be over 1!
$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('frac_problems2.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->frac_aligned_query, 0.97;

# Also, start and end should be sane
$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('frac_problems3.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->start('sbjct'), 207;
is $hit->end('sbjct'), 1051;

#WU-TBlastN test

$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('dnaEbsub_ecoli.wutblastn'));

$result = $searchio->next_result;
is($result->database_name, 'ecoli.nt');
is($result->database_letters, 4662239);
is($result->database_entries, 400);
is($result->algorithm, 'TBLASTN');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
is($result->query_name, 'gi|142865|gb|AAA22406.1|');
is($result->query_description, 'DNA primase');
is($result->query_accession, 'AAA22406.1');
is($result->query_gi, 142865);
is($result->query_length, 603);
is($result->get_parameter('matrix'), 'blosum62');

is($result->get_statistic('lambda'), '0.320');
is($result->get_statistic('kappa'), 0.136);
is($result->get_statistic('entropy'),0.387 );

is($result->get_statistic('dbentries'), 400);

@valid = ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '1.4e-73', 671]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 415);
            is($hsp->query->strand, 0);
            is($hsp->hit->start, 4778);
            is($hsp->hit->end, 6016);
            is($hsp->hit->strand, 1);
            is($hsp->length('hsp'), 421);	    
            cmp_ok($hsp->evalue,'==',1.4e-73);
            cmp_ok($hsp->pvalue,'==',1.4e-73);
            is($hsp->score,671);
            is($hsp->bits,265.8);
            is(sprintf("%.2f",$hsp->percent_identity), 35.87);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.3656);	    
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.3639);
            is(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.5297);
            is($hsp->query->frame(), 0);
            is($hsp->hit->frame(), 1);
            is($hsp->gaps('query'), 6);
            is($hsp->gaps('hit'), 8);
            is($hsp->gaps, 14);	    
            is($hsp->query_string, 'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ');
            is($hsp->hit_string, 'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ');
            is($hsp->homology_string, 'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q');
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '2 3 7-9 12 14-17 24-26 37-39 44 46 54 58 67 70-72 75-77 79-81 83 87 88 91 92 94 97 99 101-104 108 109 111-114 116 118 120 121 123 125 126 129-131 133-140 142 143 146 147 150 152 156 157 159 164-166 169 171 173-179 181-183 185 186 192 193 195 197 198 200 205 212 214 215 217 220 222 225 229 230 240 245 247 250 251 256-259 261 262 264 267 271 274-280 282 283 292 294 301 303-306 309 313 318 321 322 325 327-331 333 337-339 344 348 349 353 355 357 360 361 363 365 368 370 373-381 385-388 390-396 398 399 402 404 406-408 410');
            is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '4781-4786 4796-4804 4811-4813 4817-4828 4847-4855 4886-4894 4907-4909 4913-4915 4937-4939 4949-4951 4976-4978 4985-4993 5000-5008 5012-5020 5024-5026 5036-5041 5048-5053 5057-5059 5066-5068 5072-5077 5087-5089 5093-5101 5105-5116 5120-5122 5126-5128 5132-5137 5141-5143 5147-5152 5159-5191 5195-5200 5207-5212 5219-5221 5225-5227 5237-5242 5246-5248 5261-5269 5276-5278 5282-5284 5288-5308 5312-5320 5324-5329 5345-5350 5354-5356 5360-5368 5381-5383 5402-5404 5408-5413 5417-5419 5426-5428 5432-5434 5441-5443 5453-5458 5486-5488 5501-5503 5507-5509 5516-5521 5534-5545 5549-5554 5558-5560 5567-5569 5579-5581 5588-5608 5612-5617 5642-5644 5648-5650 5669-5671 5675-5686 5693-5695 5705-5707 5720-5722 5729-5734 5741-5743 5747-5770 5774-5776 5786-5794 5807-5809 5819-5824 5834-5836 5840-5842 5846-5848 5855-5863 5867-5869 5876-5878 5882-5884 5891-5917 5927-5938 5942-5962 5966-5971 5978-5980 5984-5986 5990-6001');
            is($hsp->ambiguous_seq_inds, 'subject');
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is($count, 1);

# WU-BLAST TBLASTX
$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('dnaEbsub_ecoli.wutblastx'));

$result = $searchio->next_result;
is($result->database_name, 'ecoli.nt');
is($result->database_letters, 4662239);
is($result->database_entries, 400);
is($result->algorithm, 'TBLASTX');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
is($result->query_name, 'gi|142864|gb|M10040.1|BACDNAE');
is($result->query_description, 'B.subtilis dnaE gene encoding DNA primase, complete cds');
is($result->query_accession, 'M10040.1');
is($result->query_gi, 142864);
is($result->query_length, 2001);
is($result->get_parameter('matrix'), 'blosum62');

is($result->get_statistic('lambda'), 0.318);
is($result->get_statistic('kappa'), 0.135);
is($result->get_statistic('entropy'),0.401 );
is($result->get_statistic('dbentries'), 400);

@valid = ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '6.4e-70', 318, 148.6],
	   [ 'gi|2367383|gb|AE000509.1|AE000509', 10589, 'AE000509', 1, 59, 29.9]
	   );
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    # using e here to deal with 0.9992 coming out right here as well
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );
    is($hit->bits, shift @$d );
    if( $count == 0 ) {
        my $hspcounter = 0;
        while( my $hsp = $hit->next_hsp ) {
            $hspcounter++;
            if( $hspcounter == 3 ) {
                # let's actually look at the 3rd HSP
                is($hsp->query->start, 441);
                is($hsp->query->end, 617);
                is($hsp->query->strand, 1);
                is($hsp->hit->start, 5192);
                is($hsp->hit->end, 5368);
                is($hsp->hit->strand, 1);
                is($hsp->length('hsp'), 59);	    
                cmp_ok($hsp->evalue,'==',6.4e-70);
                cmp_ok($hsp->pvalue,'==',6.4e-70);
                is($hsp->score,85);
                is($hsp->bits,41.8);
                is(sprintf("%.2f",$hsp->percent_identity), '32.20');
                is(sprintf("%.3f",$hsp->frac_identical('hit')), 0.322);
                is(sprintf("%.3f",$hsp->frac_identical('query')), 0.322);
                is(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.4746);
                is($hsp->query->frame(), 2);
                is($hsp->hit->frame(), 1);
                is($hsp->gaps('query'), 0);
                is($hsp->gaps('hit'), 0);
                is($hsp->gaps, 0);	    
                is($hsp->query_string, 'ALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGY');
                is($hsp->hit_string, 'ARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY');
                is($hsp->homology_string, 'A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y');
                is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '444-449 456-461 468-470 474-476 486-491 495-497 510-518 525-527 531-533 537-557 561-569 573-578 594-599 603-605 609-614');
                is(join(' ', $hsp->seq_inds('hit', 'nomatch',1)), '5195-5200 5207-5212 5219-5221 5225-5227 5237-5242 5246-5248 5261-5269 5276-5278 5282-5284 5288-5308 5312-5320 5324-5329 5345-5350 5354-5356 5360-5365');
                is($hsp->ambiguous_seq_inds, 'query/subject');
                last;
            }
        }
        is($hspcounter, 3);
    }
    elsif( $count == 1 ) {
        my $hsps_to_do = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 587);
            is($hsp->query->end, 706);
            is($hsp->query->strand, -1);
            is($hsp->hit->start, 4108);
            is($hsp->hit->end, 4227);
            is($hsp->hit->strand, -1);
            is($hsp->length('hsp'), 40);	    
            is($hsp->evalue , 7.1);
            is($hsp->pvalue , '1.00');
            is($hsp->score,59);
            is($hsp->bits,29.9);
            is(sprintf("%.2f",$hsp->percent_identity), '37.50');
            is(sprintf("%.4f",$hsp->frac_identical('hit')), '0.3750');
            is(sprintf("%.4f",$hsp->frac_identical('query')), '0.3750');
            is(sprintf("%.4f",$hsp->frac_conserved('hsp')), '0.4750');
            is($hsp->query->frame(), 2);
            is($hsp->hit->frame(), 2);
            is($hsp->gaps('query'), 0);
            is($hsp->gaps('hit'), 0);
            is($hsp->gaps, 0);
            is($hsp->query_string, 'WLPRALPEKATTAP**SWIGNMTRFLKRSKYPLPSSRLIR');
            is($hsp->hit_string, 'WLSRTTVGSSTVSPRTFWITRMKVKLSSSKVTLPSTKSTR');
            is($hsp->homology_string, 'WL R     +T +P   WI  M   L  SK  LPS++  R');
            $hsps_to_do--;
            last;
        }
        is($hsps_to_do, 0);
    }
    last if( $count++ > @valid );
}
is($count, 2);

# WU-BLAST -echofilter option test (Bug 2388)
$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('echofilter.wublastn'));

$result = $searchio->next_result;

is($result->database_name, 'NM_003201.fa');
is($result->database_letters, 1936);
is($result->database_entries, 1);
is($result->algorithm, 'BLASTN');
like($result->algorithm_version, qr/^2\.0MP\-WashU/);
like($result->query_name, qr/ref|NM_003201.1| Homo sapiens transcription factor A, mitochondrial \(TFAM\), mRNA/);
is($result->query_accession, 'NM_003201.1');

is($result->query_length, 1936);
is($result->get_statistic('lambda'), 0.192);
is($result->get_statistic('kappa'), 0.182);
is($result->get_statistic('entropy'), 0.357);
is($result->get_statistic('dbletters'), 1936);
is($result->get_statistic('dbentries'), 1);
is($result->get_parameter('matrix'), '+5,-4');

@valid = ( [ 'ref|NM_003201.1|', 1936, 'NM_003201', '0', 9680],);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 1);
            is($hsp->query->end, 1936);
            is($hsp->hit->start, 1);
            is($hsp->hit->end, 1936);
            is($hsp->length('hsp'), 1936);
            
            is($hsp->evalue , '0.');
            is($hsp->pvalue , '0.');
            is($hsp->score, 9680);
            is($hsp->bits,1458.4);	    	    
            is($hsp->percent_identity, 100);
            is($hsp->frac_identical('query'), 1.00);
            is($hsp->frac_identical('hit'), 1.00);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);


# Do a multiblast report test
$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('multi_blast.bls'));

@expected = qw(CATH_RAT CATL_HUMAN CATL_RAT PAPA_CARPA);
my $results_left = 4;
while( my $result = $searchio->next_result ) {
    is($result->query_name, shift @expected, "Multiblast query test");
    $results_left--;
}
is($results_left, 0);

# Test GCGBlast parsing

$searchio = Bio::SearchIO->new('-format' => 'blast',
			      '-file'   => test_input_file('test.gcgblast'));
$result = $searchio->next_result();

is($result->query_name, '/v0/people/staji002/test.gcg');
is($result->algorithm, 'BLASTP');
is($result->algorithm_version, '2.2.1 [Apr-13-2001]');
is($result->database_name, 'pir');
is($result->database_entries, 274514);
is($result->database_letters, 93460074);

$hit = $result->next_hit;
is($hit->description, 'F22B7.10 protein - Caenorhabditis elegans');
is($hit->name, 'PIR2:S44629');
is($hit->length, 628);
is($hit->accession, 'PIR2:S44629');
is($hit->significance, '2e-08' );
is($hit->raw_score, 136 );
is($hit->bits, '57.0' );
$hsp = $hit->next_hsp;
cmp_ok($hsp->evalue, '==', 2e-08);
is($hsp->bits, '57.0');
is($hsp->score, 136);
is(int($hsp->percent_identity), 28);
is(sprintf("%.2f",$hsp->frac_identical('query')), 0.29);
is($hsp->frac_conserved('total'), 69/135);
is($hsp->gaps('total'), 8);
is($hsp->gaps('hit'), 6);
is($hsp->gaps('query'), 2);

is($hsp->hit->start, 342);
is($hsp->hit->end, 470);
is($hsp->query->start, 3);
is($hsp->query->end, 135);

is($hsp->query_string, 'CAAEFDFMEKETPLRYTKTXXXXXXXXXXXXXXRKIISDMWGVLAKQQTHVRKHQFDHGELVYHALQLLAYTALGILIMRLKLFLTPYMCVMASLICSRQLFGW--LFCKVHPGAIVFVILAAMSIQGSANLQTQ');
is($hsp->hit_string, 'CSAEFDFIQYSTIEKLCGTLLIPLALISLVTFVFNFVKNT-NLLWRNSEEIG----ENGEILYNVVQLCCSTVMAFLIMRLKLFMTPHLCIVAALFANSKLLGGDRISKTIRVSALVGVI-AILFYRGIPNIRQQ');
is($hsp->homology_string, 'C+AEFDF++  T  +   T                 + +   +L +    +     ++GE++Y+ +QL   T +  LIMRLKLF+TP++C++A+L  + +L G   +   +   A+V VI A +  +G  N++ Q');

$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('HUMBETGLOA.tblastx'));

$result = $searchio->next_result;

isa_ok($result,'Bio::Search::Result::ResultI');
$hit = $result->next_hit;
is($hit->accession, 'AE000479');
is($hit->bits, 33.6);
$hsp = $hit->next_hsp;
is($hit->hsp->bits,$hsp->bits);

isa_ok($hsp->get_aln,'Bio::Align::AlignI');
my $writer = Bio::SearchIO::Writer::HitTableWriter->new( 
                                  -columns => [qw(
                                                  query_name
                                                  query_length
                                                  hit_name
                                                  hit_length
						  bits
						  score
                                                  frac_identical_query
                                                  expect
                                                  )]  );

my $outfile = test_output_file();
my $out = Bio::SearchIO->new(-writer => $writer,
			 -file   => ">$outfile");
$out->write_result($result, 1);
ok(-s $outfile);
my $writerhtml = Bio::SearchIO::Writer::HTMLResultWriter->new();
my $outhtml = Bio::SearchIO->new(-writer => $writerhtml,
				-file   => ">$outfile");
$outhtml->write_result($result, 1);
ok(-s $outfile);

#test all the database accession number formats
$searchio = Bio::SearchIO->new(-format => 'blast',
				 -file   => test_input_file('testdbaccnums.out'));
$result = $searchio->next_result;

@valid = ( ['pir||T14789','T14789','T14789','CAB53709','AAH01726'],
	   ['gb|NP_065733.1|CYT19', 'NP_065733','CYT19'],
	   ['emb|XP_053690.4|Cyt19','XP_053690'],
	   ['dbj|NP_056277.2|DKFZP586L0724','NP_056277'],
	   ['prf||XP_064862.2','XP_064862'],
	   ['pdb|BAB13968.1|1','BAB13968'],
	   ['sp|Q16478|GLK5_HUMAN','Q16478'],
	   ['pat|US|NP_002079.2','NP_002079'],
	   ['bbs|NP_079463.2|','NP_079463'],
	   ['gnl|db1|NP_002444.1','NP_002444'],
	   ['ref|XP_051877.1|','XP_051877'],
	   ['lcl|AAH16829.1|','AAH16829'],
	   ['gi|1|gb|NP_065733.1|CYT19','NP_065733'],
	   ['gi|2|emb|XP_053690.4|Cyt19','XP_053690'],
	   ['gi|3|dbj|NP_056277.2|DKFZP586L0724','NP_056277'],
	   ['gi|4|pir||T14789','T14789'],
	   ['gi|5|prf||XP_064862.2','XP_064862'],
	   ['gi|6|pdb|BAB13968.1|1','BAB13968'],
	   ['gi|7|sp|Q16478|GLK5_HUMAN','Q16478'],
	   ['gi|8|pat|US|NP_002079.2','NP_002079'],
	   ['gi|9|bbs|NP_079463.2|','NP_079463'],
	   ['gi|10|gnl|db1|NP_002444.1','NP_002444'],
	   ['gi|11|ref|XP_051877.1|','XP_051877'],
	   ['gi|12|lcl|AAH16829.1|','AAH16829'],
	   ['MY_test_ID','MY_test_ID']
	   );

$hit = $result->next_hit;
my $d = shift @valid;
is($hit->name, shift @$d);
is($hit->accession, shift @$d);
my @accnums = $hit->each_accession_number;
foreach my $a (@accnums) {
	is($a, shift @$d);
}
$d = shift @valid;
$hit = $result->next_hit;
is($hit->name, shift @$d);
is($hit->accession, shift @$d);
is($hit->locus, shift @$d);

$hits_left = 23;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->accession, shift @$d);
    $hits_left--;
}
is($hits_left, 0);

# Parse MEGABLAST

# parse the BLAST-like output
my $infile = test_input_file('503384.MEGABLAST.2');
my $in = Bio::SearchIO->new(-file => $infile,
			   -format => 'blast'); # this is megablast 
                                                # blast-like output
my $r = $in->next_result;
my @dcompare = ( ['Contig3700', 5631, 396, 785,'0.0', 785, '0.0', 396, 639, 12, 
		  8723,9434, 1, 4083, 4794, -1],
                 ['Contig3997', 12734, 335, 664,'0.0', 664, '0.0', 335, 401, 0, 
		  1282, 1704, 1, 1546, 1968,-1 ],
                 ['Contig634', 858, 245, 486,'1e-136', 486, '1e-136', 245, 304, 3, 
		  7620, 7941, 1, 1, 321, -1],
                 ['Contig1853', 2314, 171, 339,'1e-91',339, '1e-91', 171, 204, 0,
		  6406, 6620, 1, 1691, 1905, 1]
    );

is($r->query_name, '503384');
is($r->query_description, '11337 bp 2 contigs');
is($r->query_length, 11337);
is($r->database_name, 'cneoA.nt ');
is($r->database_letters, 17206226);
is($r->database_entries, 4935);

$hits_left = 4;
while( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->raw_score, shift @$d);
    is($hit->bits, shift @$d);
    is($hit->significance, shift @$d);
    
    my $hsp = $hit->next_hsp;
    is($hsp->bits, shift @$d);
    cmp_ok($hsp->evalue, '==', shift @$d);
    is($hsp->score, shift @$d);
    is($hsp->num_identical, shift @$d);
    is($hsp->gaps('total'), shift @$d);
    is($hsp->query->start, shift @$d);
    is($hsp->query->end, shift @$d);
    is($hsp->query->strand, shift @$d);
    is($hsp->hit->start, shift @$d);
    is($hsp->hit->end, shift @$d);
    is($hsp->hit->strand, shift @$d);
    $hits_left--;
}
is($hits_left, 0);

# parse another megablast format

$infile =  test_input_file('503384.MEGABLAST.0');

# this is megablast output type 0 
$in = Bio::SearchIO->new(-file          => $infile,
			-report_format => 0,
			-format        => 'megablast'); 
$r = $in->next_result;
@dcompare = ( 
	      ['Contig634', 7620, 7941, 1, 1, 321, -1],
	      ['Contig1853', 6406, 6620, 1, 1691, 1905, 1],  
	      ['Contig3700', 8723,9434, 1, 4083, 4794, -1],
	      ['Contig3997', 1282, 1704, 1, 1546, 1968,-1 ],
	      );

is($r->query_name, '503384');

while( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    is($hit->name, shift @$d);
    my $hsp = $hit->next_hsp;
    is($hsp->query->start, shift @$d);
    is($hsp->query->end, shift @$d);
    is($hsp->query->strand, shift @$d);
    is($hsp->hit->start, shift @$d);
    is($hsp->hit->end, shift @$d);
    is($hsp->hit->strand, shift @$d);
}
is(@dcompare, 0);


# Let's test RPS-BLAST

my $parser = Bio::SearchIO->new(-format => 'blast',
			       -file   => test_input_file('ecoli_domains.rpsblast'));

$r = $parser->next_result;
is($r->query_name, 'gi|1786183|gb|AAC73113.1|');
is($r->query_gi, 1786183);
is($r->num_hits, 7);
$hit = $r->next_hit;
is($hit->name, 'gnl|CDD|3919');
is($hit->significance, 0.064);
is($hit->bits, 28.3);
is($hit->raw_score, 63);
$hsp = $hit->next_hsp;
is($hsp->query->start, 599);
is($hsp->query->end,655);
is($hsp->hit->start,23);
is($hsp->hit->end,76);


# Test PSI-BLAST parsing

$searchio = Bio::SearchIO->new('-format' => 'blast',
			       '-file'   => test_input_file('psiblastreport.out'));

$result = $searchio->next_result;

is($result->database_name, '/home/peter/blast/data/swissprot.pr');
is($result->database_entries, 88780);
is($result->database_letters, 31984247);

is($result->algorithm, 'BLASTP');
like($result->algorithm_version, qr/^2\.0\.14/);
is($result->query_name, 'CYS1_DICDI');
is($result->query_length, 343);
is($result->get_statistic('kappa') , 0.0491);
cmp_ok($result->get_statistic('lambda'), '==', 0.270);
cmp_ok($result->get_statistic('entropy'), '==', 0.230);
is($result->get_statistic('dbletters'), 31984247);
is($result->get_statistic('dbentries'), 88780);
is($result->get_statistic('effective_hsplength'), 49);
is($result->get_statistic('effectivespace'), 8124403938);
is($result->get_parameter('matrix'), 'BLOSUM62');
is($result->get_parameter('gapopen'), 11);
is($result->get_parameter('gapext'), 1);

my @valid_hit_data = ( [ 'sp|P04988|CYS1_DICDI', 343, 'P04988', '0', 721],
		       [ 'sp|P43295|A494_ARATH', 313, 'P43295', '1e-75', 281],
		       [ 'sp|P25804|CYSP_PEA', 363, 'P25804', '1e-74', 278]);
my @valid_iter_data = ( [ 127, 127, 0, 109, 18, 0, 0, 0, 0],
			[ 157, 40, 117, 2, 38, 0, 109, 3, 5]);
my $iter_count = 0;
while( $iter = $result->next_iteration ) {
    $iter_count++;
    my $di = shift @valid_iter_data;
    is($iter->number, $iter_count);

    is($iter->num_hits, shift @$di);
    is($iter->num_hits_new, shift @$di);
    is($iter->num_hits_old, shift @$di);
    is(scalar($iter->newhits_below_threshold), shift @$di);
    is(scalar($iter->newhits_not_below_threshold), shift @$di);
    is(scalar($iter->newhits_unclassified), shift @$di);
    is(scalar($iter->oldhits_below_threshold), shift @$di);
    is(scalar($iter->oldhits_newly_below_threshold), shift @$di);
    is(scalar($iter->oldhits_not_below_threshold), shift @$di);

    my $hit_count = 0;
    if ($iter_count == 1) {
        while( $hit = $result->next_hit ) {
            my $d = shift @valid_hit_data;
            
            is($hit->name, shift @$d);
            is($hit->length, shift @$d);
            is($hit->accession, shift @$d);
            is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
            is($hit->bits, shift @$d );
            
            if( $hit_count == 1 ) {
                my $hsps_left = 1;
                while( my $hsp = $hit->next_hsp ){
                    is($hsp->query->start, 32);
                    is($hsp->query->end, 340);
                    is($hsp->hit->start, 3);
                    is($hsp->hit->end, 307);
                    is($hsp->length('hsp'), 316);
                    is($hsp->start('hit'), $hsp->hit->start);
                    is($hsp->end('query'), $hsp->query->end);
                    is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
                    cmp_ok($hsp->evalue, '==', 1e-75);
                    is($hsp->score, 712);
                    is($hsp->bits, 281);
                    is(sprintf("%.1f",$hsp->percent_identity), 46.5);
                    is(sprintf("%.4f",$hsp->frac_identical('query')), 0.4757);
                    is(sprintf("%.3f",$hsp->frac_identical('hit')), 0.482);
                    is($hsp->gaps, 18);
                    $hsps_left--;
                }
                is($hsps_left, 0);
            }
            last if( $hit_count++ > @valid_hit_data );
        }
    }
}
is(@valid_hit_data, 0);
is(@valid_iter_data, 0);

# Test filtering

$searchio = Bio::SearchIO->new( '-format' => 'blast', 
                                '-file'   => test_input_file('ecolitst.bls'),
                                '-signif' => 1e-100);

@valid = qw(gb|AAC73113.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    is($hit->name, shift @valid);
}

$searchio = Bio::SearchIO->new( '-format' => 'blast', 
                                '-file'   => test_input_file('ecolitst.bls'),
                                '-score' => 100);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1| gb|AAC76994.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    is($hit->name, shift @valid);
}
is(@valid, 0);

$searchio = Bio::SearchIO->new( '-format' => 'blast', 
                                '-file'   => test_input_file('ecolitst.bls'),
                                '-bits' => 200);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    is($hit->name, shift @valid);
}
is(@valid, 0);


my $filt_func = sub{ my $hit=shift; 
                     $hit->frac_identical('query') >= 0.31
                     };

$searchio = Bio::SearchIO->new( '-format' => 'blast', 
                                '-file'   => test_input_file('ecolitst.bls'),
                                '-hit_filter' => $filt_func);

@valid = qw(gb|AAC73113.1| gb|AAC76994.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    is($hit->name, shift @valid);
}
is(@valid, 0);




# bl2seq parsing testing

# this is blastp bl2seq
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.out'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, '');
is($result->algorithm, 'BLASTP');
$hit = $result->next_hit;
is($hit->name, 'ALEU_HORVU');
is($hit->length, 362);
$hsp = $hit->next_hsp;
is($hsp->score, 481);
is($hsp->bits, 191);
is(int $hsp->percent_identity, 34);
cmp_ok($hsp->evalue, '==', 2e-53);
is(int($hsp->frac_conserved*$hsp->length), 167);
is($hsp->query->start, 28);
is($hsp->query->end, 343);
is($hsp->hit->start, 60);
is($hsp->hit->end,360);
is($hsp->gaps, 27);

# this is blastn bl2seq 
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.blastn.rev'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, '');
is($result->algorithm, 'BLASTN');
is($result->query_length, 180);
$hit = $result->next_hit;
is($hit->length, 179);
is($hit->name, 'human');
$hsp = $hit->next_hsp;
is($hsp->score, 27);
is($hsp->bits, '54.0');
is(int $hsp->percent_identity, 88);
cmp_ok($hsp->evalue, '==', 2e-12);
is(int($hsp->frac_conserved*$hsp->length), 83);
is($hsp->query->start, 94);
is($hsp->query->end, 180);
is($hsp->query->strand, 1);
is($hsp->hit->strand, -1);
is($hsp->hit->start, 1);
is($hsp->hit->end,94);
is($hsp->gaps, 7);

# this is blastn bl2seq 
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.blastn'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, '');
is($result->query_length, 180);
is($result->algorithm, 'BLASTN');
$hit = $result->next_hit;
is($hit->name, 'human');
is($hit->length, 179);
$hsp = $hit->next_hsp;
is($hsp->score, 27);
is($hsp->bits, '54.0');
is(int $hsp->percent_identity, 88);
cmp_ok($hsp->evalue,'==', 2e-12);
is(int($hsp->frac_conserved*$hsp->length), 83);
is($hsp->query->start, 94);
is($hsp->query->end, 180);
is($hsp->query->strand, 1);
is($hsp->hit->strand, 1);
is($hsp->hit->start, 86);
is($hsp->hit->end,179);
is($hsp->gaps, 7);

# this is blastp bl2seq
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.bug940.out'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, 'zinc');
is($result->algorithm, 'BLASTP');
is($result->query_description, 'finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0');
is($result->query_length, 469);
$hit = $result->next_hit;
is($hit->name, 'gi|4507985|');
is($hit->ncbi_gi, 4507985);
is($hit->description,'zinc finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0');
is($hit->length, 469);
$hsp = $hit->next_hsp;
is($hsp->score, 1626);
is($hsp->bits, 637);
is(int $hsp->percent_identity, 66);
cmp_ok($hsp->evalue, '==', 0.0);
is(int($hsp->frac_conserved*$hsp->length), 330);
is($hsp->query->start, 121);
is($hsp->query->end, 469);
is($hsp->hit->start, 1);
is($hsp->hit->end,469);
is($hsp->gaps, 120);
ok($hit->next_hsp); # there is more than one HSP here, 
                    # make sure it is parsed at least

# cannot distinguish between blastx and tblastn reports
# so we're only testing a blastx report for now

# this is blastx bl2seq
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.blastx.out'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, 'AE000111.1');
is($result->query_description, 'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome');
is($result->algorithm, 'BLASTX');
is($result->query_length, 720);
$hit = $result->next_hit;
is($hit->name, 'AK1H_ECOLI');
is($hit->description,'P00561 Bifunctional aspartokinase/homoserine dehydrogenase I (AKI-HDI) [Includes: Aspartokinase I ; Homoserine dehydrogenase I ]');
is($hit->length, 820);
$hsp = $hit->next_hsp;
is($hsp->score, 634);
is($hsp->bits, 248);
is(int $hsp->percent_identity, 100);
cmp_ok($hsp->evalue, '==' ,2e-70);
is(int($hsp->frac_conserved*$hsp->length), 128);
is($hsp->query->start, 1);
is($hsp->query->end, 384);
is($hsp->hit->start, 1);
is($hsp->hit->end,128);
is($hsp->gaps, 0);
is($hsp->query->frame,0);
is($hsp->hit->frame,0);
is($hsp->query->strand,-1);
is($hsp->hit->strand,0);

# this is tblastx bl2seq (self against self)
$searchio = Bio::SearchIO->new(-format => 'blast',
			      -file   => test_input_file('bl2seq.tblastx.out'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->query_name, 'Escherichia');
is($result->algorithm, 'TBLASTX');
is($result->query_description, 'coli K-12 MG1655 section 1 of 400 of the complete genome');
is($result->query_length, 720);
$hit = $result->next_hit;
is($hit->name, 'gi|1786181|gb|AE000111.1|AE000111');
is($hit->ncbi_gi, 1786181);
is($hit->description,'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome');
is($hit->length, 720);
$hsp = $hit->next_hsp;
is($hsp->score, 1118);
is($hsp->bits, 515);
is(int $hsp->percent_identity, 95);
cmp_ok($hsp->evalue, '==', 1e-151);
is(int($hsp->frac_conserved*$hsp->length), 229);
is($hsp->query->start, 1);
is($hsp->query->end, 720);
is($hsp->hit->start, 1);
is($hsp->hit->end,720);
is($hsp->gaps, 0);
is($hsp->query->frame,0);
is($hsp->hit->frame,0);
is($hsp->query->strand,1);
is($hsp->hit->strand,1);

# this is NCBI tblastn
$searchio = Bio::SearchIO->new(-format => 'blast',
										-file   => test_input_file('tblastn.out'));
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
is($result->algorithm, 'TBLASTN');
$hit = $result->next_hit;
is($hit->name,'gi|10040111|emb|AL390796.6|AL390796');

# test blasttable output
my @eqset = qw( c200-vs-yeast.BLASTN.m9);
$searchio = Bio::SearchIO->new(-file => test_input_file('c200-vs-yeast.BLASTN'),
			      -format => 'blast');
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
my %ref = &result2hash($result);
is( scalar keys %ref, 67);
$searchio = Bio::SearchIO->new(-file => test_input_file('c200-vs-yeast.BLASTN.m8'),
			      -program_name => 'BLASTN',
			      -format => 'blasttable');
$result = $searchio->next_result;
my %tester = &result2hash($result);
is( scalar keys %tester, 67);
foreach my $key ( sort keys %ref ) {
    is($tester{$key}, $ref{$key},$key);
}      

# test WU-BLAST blasttable output
$searchio = Bio::SearchIO->new(-file => test_input_file('test1.wublastp'),
			      			   -format => 'blast');
$result = $searchio->next_result;
isa_ok($result,'Bio::Search::Result::ResultI');
my %wuref = &result2hash($result);
is( scalar keys %wuref, 31);
$searchio = Bio::SearchIO->new(-file => test_input_file('test1.blasttab3'),
			      -program_name => 'BLASTP',
			      -format => 'blasttable');
$result = $searchio->next_result;
my %wutester = &result2hash($result);
is( scalar keys %wutester, 31);
foreach my $key ( sort keys %ref ) {
    is($wutester{$key}, $wuref{$key},$key);
}      

# Test Blast parsing with B=0 (WU-BLAST)
$searchio = Bio::SearchIO->new(-file   => test_input_file('no_hsps.blastp'),
			      -format => 'blast');
$result = $searchio->next_result;
is($result->query_name, 'mgri:MG00189.3');
$hit = $result->next_hit;
is($hit->name, 'mgri:MG00189.3');
is($hit->description, 'hypothetical protein 6892 8867 +');
is($hit->bits, 3098);
is($hit->significance, '0.');

$hit = $result->next_hit;
is($hit->name, 'fgram:FG01141.1');
is($hit->description, 'hypothetical protein 47007 48803 -');
is($hit->bits, 2182);
is($hit->significance, '4.2e-226');
is($result->num_hits, 415);
# Let's now test if _guess_format is doing its job correctly
my %pair = ( 'filename.blast'  => 'blast',
	     'filename.bls'    => 'blast',
	     'f.blx'           => 'blast',
	     'f.tblx'          => 'blast',
	     'fast.bls'        => 'blast',
	     'f.fasta'         => 'fasta',
	     'f.fa'            => 'fasta',
	     'f.fx'            => 'fasta',
	     'f.fy'            => 'fasta',
	     'f.ssearch'       => 'fasta',
	     'f.SSEARCH.m9'    => 'fasta',
	     'f.m9'            => 'fasta',
	     'f.psearch'       => 'fasta',
	     'f.osearch'       => 'fasta',
	     'f.exon'          => 'exonerate',
	     'f.exonerate'     => 'exonerate',
	     'f.blastxml'      => 'blastxml',
	     'f.xml'           => 'blastxml');
while( my ($file,$expformat) = each %pair ) {
    is(Bio::SearchIO->_guess_format($file),$expformat, "$expformat for $file");
}


# Test Wes Barris's reported bug when parsing blastcl3 output which
# has integer overflow

$searchio = Bio::SearchIO->new(-file => test_input_file('hsinsulin.blastcl3.blastn'),
			      -format => 'blast');
$result = $searchio->next_result;
is($result->query_name, 'human');
is($result->database_letters(), '-24016349'); 
# this is of course not the right length, but is the what blastcl3 
# reports, the correct value is
is($result->get_statistic('dbletters'),'192913178');
is($result->get_statistic('dbentries'),'1867771');


# test for links and groups being parsed out of WU-BLAST properly
$searchio = Bio::SearchIO->new(-format => 'blast',
			       -file   => test_input_file('brassica_ATH.WUBLASTN'));
ok($result = $searchio->next_result);
ok($hit = $result->next_hit);
ok($hsp = $hit->next_hsp);
is($hsp->links,'(1)-3-2');
is($hsp->query->strand, 1);
is($hsp->hit->strand, 1);
is($hsp->hsp_group, '1');

## Web blast result parsing

$searchio = Bio::SearchIO->new(-format => 'blast',
			       -file   => test_input_file('catalase-webblast.BLASTP'));
ok($result = $searchio->next_result);
ok($hit = $result->next_hit);
is($hit->name, 'gi|40747822|gb|EAA66978.1|', 'full hit name');
is($hit->accession, 'EAA66978', 'hit accession');
is($hit->ncbi_gi, 40747822);
ok($hsp = $hit->next_hsp);
is($hsp->query->start, 1, 'query start');
is($hsp->query->end, 528, 'query start');

# tests for new BLAST 2.2.13 output
$searchio = Bio::SearchIO->new(-format => 'blast',
							  -file   => test_input_file('new_blastn.txt'));

$result = $searchio->next_result;
is($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)');
is($result->database_entries, 3742891);
is($result->database_letters, 16670205594);
is($result->algorithm, 'BLASTN');
is($result->algorithm_version, '2.2.13 [Nov-27-2005]');
is($result->query_name, 'pyrR,');
is($result->query_length, 558);
is($result->get_statistic('kappa'), '0.711');
is($result->get_statistic('kappa_gapped'), '0.711');
is($result->get_statistic('lambda'), '1.37');
is($result->get_statistic('lambda_gapped'), '1.37');
is($result->get_statistic('entropy'), '1.31');
is($result->get_statistic('entropy_gapped'), '1.31');
is($result->get_statistic('dbletters'), '-509663586');
is($result->get_statistic('dbentries'), 3742891);
is($result->get_statistic('effective_hsplength'), undef);
is($result->get_statistic('effectivespace'), undef);
is($result->get_parameter('matrix'), 'blastn matrix:1 -3');
is($result->get_parameter('gapopen'), 5);
is($result->get_parameter('gapext'), 2);
is($result->get_statistic('S2'), '60');
is($result->get_statistic('S2_bits'), '119.4');
is($result->get_parameter('expect'), '1e-23');
is($result->get_statistic('num_extensions'), '117843');


@valid = ( [ 'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', 41400296, '6e-059', 119, 236],
	      [ 'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', 54013472, '4e-026', 64, 127],
	      [ 'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', 57546753, '1e-023', 60, 119]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->ncbi_gi, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    is($hit->raw_score, shift @$d );
    is($hit->bits, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            is($hsp->query->start, 262);
            is($hsp->query->end, 552);
            is($hsp->hit->start, 1166897);
            is($hsp->hit->end, 1167187);
            is($hsp->length('hsp'), 291);
            is($hsp->start('hit'), $hsp->hit->start);
            is($hsp->end('query'), $hsp->query->end);
            is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            cmp_ok($hsp->evalue, '==', 6e-59);
            is($hsp->score, 119);
            is($hsp->bits,236);	    	    
            is(sprintf("%.2f",$hsp->percent_identity), 85.22);
            is(sprintf("%.4f",$hsp->frac_identical('query')), 0.8522);
            is(sprintf("%.4f",$hsp->frac_identical('hit')), 0.8522);
            is($hsp->gaps, 0);
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);

# Bug 2189
$searchio = Bio::SearchIO->new(-format => 'blast',
							  -file   => test_input_file('blastp2215.blast'));

$result = $searchio->next_result;
is($result->database_entries, 4460989);
is($result->database_letters, 1533424333);
is($result->algorithm, 'BLASTP');
is($result->algorithm_version, '2.2.15 [Oct-15-2006]');
is($result->query_name, 'gi|15608519|ref|NP_215895.1|');
is($result->query_gi, 15608519);
is($result->query_length, 193);
@hits = $result->hits;
is(scalar(@hits), 10);
is($hits[1]->accession,'1W30');
is($hits[4]->significance,'2e-72');
is($hits[7]->bits,'254');
$result = $searchio->next_result;
is($result->database_entries, 4460989);
is($result->database_letters, 1533424333);
is($result->algorithm, 'BLASTP');
is($result->algorithm_version, '2.2.15 [Oct-15-2006]');
is($result->query_name, 'gi|15595598|ref|NP_249092.1|');
is($result->query_length, 423);
@hits = $result->hits;
is(scalar(@hits), 10);
is($hits[1]->accession,'ZP_00972546');
is($hits[2]->ncbi_gi, 116054132);
is($hits[4]->significance, '0.0');
is($hits[7]->bits, 624);

# Bug 2246
$searchio = Bio::SearchIO->new(-format => 'blast',
                               -verbose => -1,
							  -file   => test_input_file('bug2246.blast'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->name, 'UniRef50_Q9X0H5';
is $hit->length, 0;
is $hit->accession, 'UniRef50_Q9X0H5';
is $hit->description, 'Cluster: Histidyl-tRNA synthetase; n=4; Thermoto...';
is $hit->bits, 23;
is $hit->significance, 650;

# Bug 1986
$searchio = Bio::SearchIO->new(-format => 'blast',
                               -verbose => -1,
							  -file   => test_input_file('bug1986.blastp'));
$result = $searchio->next_result;
$hit = $result->next_hit;
is $hit->name, 'ENSP00000350182';
is $hit->length, 425;
is $hit->accession, 'ENSP00000350182';
is $hit->description, 'pep:novel clone::BX322644.8:4905:15090:-1 gene:ENSG00000137397 transcript:ENST00000357569 ';
is $hit->raw_score, 301;
is $hit->bits, 120;
is $hit->significance, 3e-27;
$hit = $result->next_hit;
is $hit->name, 'ENSP00000327738';
is $hit->length, 468;
is $hit->accession, 'ENSP00000327738';
is $hit->description, 'pep:known-ccds chromosome:NCBI36:4:189297592:189305643:1 gene:ENSG00000184108 transcript:ENST00000332517 CCDS3851.1';
is $hit->raw_score, 289;
is $hit->bits, 115;
is $hit->significance, 8e-26;

# Bug 1986, pt. 2

# handle at least the first iteration with BLAST searches using databases
# containing non-unique IDs

my $file = test_input_file('bug1986.blast2');
my %unique_accs;
open (my $IN,$file) or die $!;

while (<$IN>) {
  last if (/^Sequences/);
}
$count = 1;
while (<$IN>) {
  chomp;
  next if m{^\s*$};
  next unless ($_);
  last if m{^>};
  my ($accession) = split(/\s+/);
  #print "Real Hit $count = $accession\n";
  $unique_accs{$accession}++;
  #last if ($count == 10);
  ++$count;
}
close ($IN);

is ($count, 495);
is (scalar(keys %unique_accs), 490);

my %search_accs;

$searchio = Bio::SearchIO->new(-format => 'blast',
                               -verbose => -1,
							  -file   => $file);
$result = $searchio->next_result;
$count = 1;
while (my $hit = $result->next_hit) {
	$search_accs{$hit->accession}++;
	$count++;
}

is ($count, 495);
is (scalar(keys %search_accs), 490);

is_deeply(\%unique_accs, \%search_accs);

# bug 2391 - long query names

$file = test_input_file('bug2391.megablast');

$searchio = Bio::SearchIO->new(-format => 'blast',
							  -file   => $file);
$result = $searchio->next_result;

# data is getting munged up with long names
is($result->query_name, 'c6_COX;c6_QBL;6|31508172;31503325;31478402|rs36223351|1|dbSNP|C/G');
is($result->query_description, '');

# bug 2399 - catching Expect(n) values

$file = test_input_file('bug2399.tblastn');

$searchio = Bio::SearchIO->new(-format => 'blast',
							  -file   => $file);
my $total_n = 0;
while(my $query = $searchio->next_result) {
    while(my $subject = $query->next_hit) {
        $total_n += grep{$_->n} $subject->hsps;
    }
}
is($total_n, 10);

# bug 2473 - fasta3.4 parsing with -U option

$file = test_input_file('bug2473.fasta');

$searchio = Bio::SearchIO->new(-format => 'fasta',
							  -file   => $file);

while(my $res = $searchio->next_result) {
    is($res->query_name, 'total:39860_L:12096_-3:12346_0:617_+3:14801');
    is($res->query_description, '');
    is($res->query_length, 22);
    is($res->algorithm, 'FASTN');
}

# BLAST 2.2.18+ tabular output (has 13 columns instead of 12)
$file = test_input_file('2008.blasttable');

$searchio = Bio::SearchIO->new(-format => 'blasttable',
							  -file   => $file);

while(my $res = $searchio->next_result) {
    is($res->query_name, 'gi|1786183|gb|AAC73113.1|');
    is($res->algorithm, 'BLASTP');
    is($res->algorithm_version, '2.2.18+');
    my $hit = $res->next_hit;
    is($hit->name, 'gi|34395933|sp|P00561.2|AK1H_ECOLI');
    $hit = $res->next_hit;
    my $hsp = $hit->next_hsp;
    is($hsp->bits, 331);
    is($hsp->evalue, '2e-91');
    is($hsp->start('hit'), 16);
    is($hsp->end('hit'), 805);
    is($hsp->start('query'), 5);
    is($hsp->end('query'), 812);
    is($hsp->length, 821);
    is($hsp->gaps, 14);
}


# some utilities
# a utility function for comparing result objects
sub result2hash {
    my ($result) = @_;
    my %hash;
    $hash{'query_name'} = $result->query_name;
    my $hitcount = 1;
    my $hspcount = 1;
    foreach my $hit ( $result->hits ) {
	$hash{"hit$hitcount\_name"}   =  $hit->name;
	# only going to test order of magnitude
	# too hard as these don't always match
#	$hash{"hit$hitcount\_signif"} =  
#	    ( sprintf("%.0e",$hit->significance) =~ /e\-?(\d+)/ );
	$hash{"hit$hitcount\_bits"}   =  sprintf("%d",$hit->bits);

	foreach my $hsp ( $hit->hsps ) {
	    $hash{"hsp$hspcount\_bits"}   = sprintf("%d",$hsp->bits);
	    # only going to test order of magnitude
 	    # too hard as these don't always match
#	    $hash{"hsp$hspcount\_evalue"} =  
#		( sprintf("%.0e",$hsp->evalue) =~ /e\-?(\d+)/ );
	    $hash{"hsp$hspcount\_qs"}     = $hsp->query->start;
	    $hash{"hsp$hspcount\_qe"}     = $hsp->query->end;
	    $hash{"hsp$hspcount\_qstr"}   = $hsp->query->strand;
	    $hash{"hsp$hspcount\_hs"}     = $hsp->hit->start;
	    $hash{"hsp$hspcount\_he"}     = $hsp->hit->end;
	    $hash{"hsp$hspcount\_hstr"}   = $hsp->hit->strand;

	    #$hash{"hsp$hspcount\_pid"}     = sprintf("%d",$hsp->percent_identity);
	    #$hash{"hsp$hspcount\_fid"}     = sprintf("%.2f",$hsp->frac_identical);
	    $hash{"hsp$hspcount\_gaps"}    = $hsp->gaps('total');
	    $hspcount++;
	}
	$hitcount++;
    }
    return %hash;
}

__END__

Useful for debugging:

    if ($iter_count == 3) {
	print "NEWHITS:\n";
	foreach ($iter->newhits) {
	    print "  " . $_->name . "\n";
	}
	print "\nOLDHITS:\n";
	foreach ($iter->oldhits) {
	    print "  " . $_->name . "\n";
	}
	print "\nNEWHITS BELOW:\n";
	foreach ($iter->newhits_below_threshold) {
	    print "  " . $_->name . "\n";
	}
	print "\nNEWHITS NOT BELOW:\n";
	foreach ($iter->newhits_not_below_threshold) {
	    print "  " . $_->name . "\n";
	}
	print "\nNEWHITS UNCLASSIFIED:\n";
	foreach ($iter->newhits_unclassified) {
	    print "  " . $_->name . "\n";
	}
	print "\nOLDHITS BELOW:\n";
	foreach ($iter->oldhits_below_threshold) {
	    print "  " . $_->name . "\n";
	}
	print "\nOLDHITS NEWLY BELOW:\n";
	foreach ($iter->oldhits_newly_below_threshold) {
	    print "  " . $_->name . "\n";
	}
	print "\nOLDHITS NOT BELOW:\n";
	foreach ($iter->oldhits_not_below_threshold) {
	    print "  " . $_->name . "\n";
	}
    }
