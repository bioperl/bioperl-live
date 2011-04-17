# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 391,
			   -requires_module => 'XML::SAX');
	
	use_ok('Bio::SearchIO');
}

my ($searchio, $result,$iter,$hit,$hsp);

# XML encoding/decoding done within XML::SAX now, though some parsers
# do not work properly (XML::SAX::PurePerl, XML::LibXML::SAX)

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

SKIP: {
	# this should be fixed with newer installations of XML::SAX::Expat, but as we
	# don't require a certain version (multiple backends can be used) we catch
	# and skip if needed 
	if ($@ && $@ =~ m{Handler could not resolve external entity}) {
		skip("Older versions of XML::SAX::Expat may not work with XML tests; skipping",297);
	} elsif ($@) {
		skip("Problem with XML::SAX setup: $@. Check ParserDetails.ini; skipping XML tests",297);
	}
	is($searchio->result_count, 1);
	
	# basic ResultI data
	isa_ok($result, 'Bio::Search::Result::ResultI');
	is($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam', 'database_name()');
	is($result->query_name,'gi|1786182|gb|AAC73112.1|','query_name()');
	is($result->query_description, '(AE000111) thr operon leader peptide [Escherichia coli]');
	is($result->query_accession, 'AAC73112.1');
	is($result->query_gi, 1786182);
	is($result->query_length, 21);
	is($result->algorithm, 'BLASTP');
	is($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');
	
	# check parameters
	is($result->available_parameters, 8);
	is($result->get_parameter('matrix'), 'BLOSUM62');
	float_is($result->get_parameter('expect'), '1e-05');
	is($result->get_parameter('include'), 0);
	is($result->get_parameter('match'), 0);
	is($result->get_parameter('mismatch'), 0);
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), 'F');
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 0);
	is($result->database_letters, 0);
	is($result->get_statistic('hsplength'), 0);
	float_is($result->get_statistic('effectivespace'), 4.16497e+11);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.041);
	is($result->get_statistic('entropy'), 0.14);
	
	# this result actually has a hit
	$result = $searchio->next_result;
	
	# does the parser catch everything in the next result?
	is($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam', 'database_name()');
	is($result->query_name,'gi|1786183|gb|AAC73113.1|');
	is($result->query_description, '(AE000111) aspartokinase I, homoserine dehydrogenase I [Escherichia coli]');
	is($result->query_accession, 'AAC73113.1');
	is($result->query_gi, 1786183);
	is($result->query_length, 820);
	is($result->algorithm, 'BLASTP');
	is($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');	
	
	is($searchio->result_count, 2);

	# check parameters
	is($result->available_parameters, 8);
	is($result->get_parameter('matrix'), 'BLOSUM62');
	float_is($result->get_parameter('expect'), '1e-05');
	is($result->get_parameter('include'), 0);
	is($result->get_parameter('match'), 0);
	is($result->get_parameter('mismatch'), 0);
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), 'F');
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 0);
	is($result->database_letters, 0);
	is($result->get_statistic('hsplength'), 0);
	float_is($result->get_statistic('effectivespace'), 3.82682e+07);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.041);
	is($result->get_statistic('entropy'), 0.14);
	
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
	float_is($hsp->evalue, 1.46134e-90);
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
	# gaps should match calculated sequence indices for gaps and vice versa
	is($hsp->gaps('total'), $hsp->seq_inds('hit', 'gaps') + $hsp->seq_inds('query', 'gaps'));
	is($hsp->gaps('hit'), $hsp->seq_inds('hit', 'gaps'));
	is($hsp->gaps('query'), $hsp->seq_inds('query', 'gaps'));
	is($hsp->length('total'), 326);
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
	is($searchio->result_count, 1);
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
	is($searchio->result_count, 1);
	is($result->database_name,'Hs15_up1000');
	is($result->query_name,'NM_011441_up_1000_chr1_4505586_r');
	is($result->query_description,'chr1:4505586-4506585');
	is($result->query_accession,'NM_011441_up_1000_chr1_4505586_r');
	is($result->query_gi, '');
	is($result->query_length,'1000');
	
	# check parameters
	is($result->available_parameters, 6);
	is($result->get_parameter('matrix'), undef); # not set
	float_is($result->get_parameter('expect'), 10);
	is($result->get_parameter('include'), undef); # not set
	is($result->get_parameter('match'), 1);
	is($result->get_parameter('mismatch'), -3);
	is($result->get_parameter('gapopen'), 5);
	is($result->get_parameter('gapext'), 2);
	is($result->get_parameter('filter'), 'D');
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 17516);
	is($result->database_letters, 17516000);
	is($result->get_statistic('hsplength'), 0);
	float_is($result->get_statistic('effectivespace'), 1.69255e+10);
	is($result->get_statistic('lambda'), 1.37407);
	is($result->get_statistic('kappa'), 0.710605);
	is($result->get_statistic('entropy'), 1.30725);
	
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
	is($searchio->result_count, 1);
	is($result->database_name,'nr');
	is($result->algorithm,'BLASTP');
	is($result->algorithm_version,'BLASTP 2.2.15 [Oct-15-2006]');
	# some XML::SAX parsers (PurePerl, XML::SAX::LibXML) don't decode entities
	# properly, not fixable using decode_entities()
	like($result->algorithm_reference, qr{Nucleic Acids Res} ); 
	is($result->query_name,'gi|15600734|ref|NP_254228.1|');
	is($result->query_description,'dihydroorotase [Pseudomonas aeruginosa PAO1]');
	is($result->query_accession,'NP_254228.1');
	is($result->query_gi, 15600734);
	is($result->query_length,'445');	

	# check parameters
	is($result->available_parameters, 4);
	is($result->get_parameter('matrix'), 'BLOSUM62'); 
	float_is($result->get_parameter('expect'), 10);
	is($result->get_parameter('include'), undef); # not set
	is($result->get_parameter('match'), undef);   # not set
	is($result->get_parameter('mismatch'), undef);# not set
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), undef);  # not set
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 4299737);
	is($result->database_letters, 1479795817);
	is($result->get_statistic('hsplength'), 0);
	float_is($result->get_statistic('effectivespace'), 0);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.041);
	is($result->get_statistic('entropy'), 0.14);
	
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
	float_is($hsp->evalue, 0);
	is($hsp->score, 2251);
	is($hsp->bits,871.692);
	is($hsp->query->start, 1);
	is($hsp->query->end,445);
	is($hsp->hit->start, 1);
	is($hsp->hit->end, 445);
	is($hsp->query->frame,0);
	is($hsp->hit->frame,0);
	
	$result = $searchio->next_result;
	is($searchio->result_count, 2);
	is($result->database_name,'nr'); 
	is($result->algorithm,'BLASTP');
	is($result->algorithm_version,'BLASTP 2.2.15 [Oct-15-2006]'); 
	like($result->algorithm_reference, qr{Nucleic Acids Res} );
	is($result->query_name,'gi|15598723|ref|NP_252217.1|');
	is($result->query_description,'dihydroorotase [Pseudomonas aeruginosa PAO1]');
	is($result->query_accession,'NP_252217.1');
	is($result->query_gi, 15598723);
	is($result->query_length,'348');
	
	# check parameters
	is($result->available_parameters, 4);
	is($result->get_parameter('matrix'), 'BLOSUM62'); 
	float_is($result->get_parameter('expect'), 10);
	is($result->get_parameter('include'), undef); # not set
	is($result->get_parameter('match'), undef);   # not set
	is($result->get_parameter('mismatch'), undef);# not set
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), undef);  # not set
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 4299737);
	is($result->database_letters, 1479795817);
	is($result->get_statistic('hsplength'), 0);
	float_is($result->get_statistic('effectivespace'), 0);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.041);
	is($result->get_statistic('entropy'), 0.14);

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
	float_is($hsp->evalue, 0);
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
	
	$result = $searchio->next_result;
	is($searchio->result_count, 1);    
	is($result->database_name, 'AL591824.faa');
	is($result->algorithm, 'BLASTP');
	like($result->algorithm_version, qr/2\.2\.16/);
	is($result->query_name, 'gi|1373160|gb|AAB57770.1|');
	is($result->query_accession, 'AAB57770.1');
	is($result->query_gi, '1373160');
	is($result->query_length, 173);
	
	# check parameters
	is($result->available_parameters, 6);
	is($result->get_parameter('matrix'), 'BLOSUM62'); 
	float_is($result->get_parameter('expect'), 10);
	is($result->get_parameter('include'), 0.002); 
	is($result->get_parameter('match'), undef);   # not set
	is($result->get_parameter('mismatch'), undef);# not set
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), 'F');  
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 2846);
	is($result->database_letters, 870878);
	is($result->get_statistic('hsplength'), 75);
	float_is($result->get_statistic('effectivespace'), 6.44279e+07);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.0475563);
	is($result->get_statistic('entropy'), 0.14);
	
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
				float_is($hit->significance, shift @$d);
				is($hit->bits, shift @$d );
				if( $hit_count == 1 ) {
					my $hsps_left = 1;
					while( my $hsp = $hit->next_hsp ){
						is($hsp->query->start, 4);
						is($hsp->query->end, 155);
						is($hsp->hit->start, 475);
						is($hsp->hit->end, 617);
						is($hsp->length('total'), 153);
						is($hsp->start('hit'), $hsp->hit->start);
						is($hsp->end('query'), $hsp->query->end);
						is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
						float_is($hsp->evalue, 0.000286309);
						is($hsp->score, 86);
						is($hsp->bits, 37.7354);
						is(sprintf("%.1f",$hsp->percent_identity), 20.9);
						is(sprintf("%.4f",$hsp->frac_identical('query')), 0.2105);
						is(sprintf("%.3f",$hsp->frac_identical('hit')), 0.224);
						is($hsp->gaps('total'), 11);
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
	is($searchio->result_count, 2);    
	is($result->database_name, 'AL591824.faa');
	is($result->algorithm, 'BLASTP');
	like($result->algorithm_version, qr/2\.2\.16/);
	is($result->query_name, 'gi|154350371|gb|ABS72450.1|');
	is($result->query_accession, 'ABS72450.1');
	is($result->query_gi, '154350371');
	is($result->query_length, 378);

	# check parameters
	is($result->available_parameters, 6);
	is($result->get_parameter('matrix'), 'BLOSUM62'); 
	float_is($result->get_parameter('expect'), 10);
	is($result->get_parameter('include'), 0.002); 
	is($result->get_parameter('match'), undef);   # not set
	is($result->get_parameter('mismatch'), undef);# not set
	is($result->get_parameter('gapopen'), 11);
	is($result->get_parameter('gapext'), 1);
	is($result->get_parameter('filter'), 'F');  
	
	# check statistics
	is($result->available_statistics, 5);
	is($result->database_entries, 2846);
	is($result->database_letters, 870878);
	is($result->get_statistic('hsplength'), 82);
	float_is($result->get_statistic('effectivespace'), 1.88702e+08);
	is($result->get_statistic('lambda'), 0.267);
	is($result->get_statistic('kappa'), 0.0450367);
	is($result->get_statistic('entropy'), 0.14);
	
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
				float_is($hit->significance, shift @$d);
				is($hit->bits, shift @$d );
				if( $hit_count == 1 ) {
					my $hsps_left = 1;
					while( my $hsp = $hit->next_hsp ){
						is($hsp->query->start, 63);
						is($hsp->query->end, 181);
						is($hsp->hit->start, 304);
						is($hsp->hit->end, 432);
						is($hsp->length('total'), 129);
						is($hsp->start('hit'), $hsp->hit->start);
						is($hsp->end('query'), $hsp->query->end);
						is($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
						float_is($hsp->evalue, 0.0242028);
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
