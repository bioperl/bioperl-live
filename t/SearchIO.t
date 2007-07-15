# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 1443);
	
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
			   '-verbose' => -1);
		# PurePerl works with these BLAST reports, so removed verbose promotion
		$result = $searchio->next_result;
        die if !defined $result;
	};
	if ($@ && $@ =~ m{Handler couldn't resolve external entity}) {
		skip("XML::SAX::Expat does not work with XML tests; skipping",129);
	} elsif ($@) {
		skip("Problem with XML::SAX setup: $@. Check ParserDetails.ini; skipping XML tests",129);
	}
	
    isa_ok($result, 'Bio::Search::Result::ResultI');
    is($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam', 'database_name()');
    is($result->query_name,'gi|1786182|gb|AAC73112.1|','query_name()');
    is($result->query_description, '(AE000111) thr operon leader peptide [Escherichia coli]');
    is($result->query_accession, 'AAC73112.1');
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
    is(join(' ', $hsp->seq_inds('query', 'gap',1)), '532 548-551 562 649 690');
	
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
    is($result->query_length,'1000');
    $hit = $result->next_hit;
    is($hit->name,'NM_001938_up_1000_chr1_93161154_f');
    is($hit->description,'chr1:93161154-93162153');
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
    is($result->query_length,'348');
    $hit = $result->next_hit;
    is($hit->name,'gi|15598723|ref|NP_252217.1|');
    is($hit->description,'dihydroorotase [Pseudomonas aeruginosa PAO1] '.
       '>gi|6226683|sp|P72170|PYRC_PSEAE Dihydroorotase (DHOase) '.
       '>gi|9949676|gb|AAG06915.1|AE004773_4 dihydroorotase [Pseudomonas aeruginosa PAO1] '.
       '>gi|3868712|gb|AAC73109.1| dihydroorotase [Pseudomonas aeruginosa]');
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


my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 1567],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922', '1e-91', 332],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994', '3e-47', 184]);
my $count = 0;
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

@valid = ( [ 'gb|AE000479.1|AE000479', 10934, 'AE000479', '0.13', 34],
	   [ 'gb|AE000302.1|AE000302', 10264, 'AE000302', '0.61', 31],
	   [ 'gb|AE000277.1|AE000277', 11653, 'AE000277', '0.84', 31]);
$count = 0;

while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is($hit->significance, shift @$d );
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
            is(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '355 364 365 367 368 370 371 373-375');
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


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309', 1.2, 29.2],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148', 2.1, 27.4],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494',  2.1, 25.9]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
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
            $hsps_left--;
        }
        is($hsps_left, 0);
    }
    last if( $count++ > @valid );
}
is(@valid, 0);
is($result->hits, 58);
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

@valid = ( [ 'gb|AY052359.1|', 2826, 'AY052359', '3e-18', 96, 1, 60, 
	     '1.0000'],
	   [ 'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 96, 1, 60, 
	     '1.0000' ],
	   [ 'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.04', 42, 35, 55, 
	     '0.3500']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    is($hit->name, shift @$d);
    is($hit->length, shift @$d);
    is($hit->accession, shift @$d);
    is(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
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
            is(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity), 96.67);
            is(sprintf("%.2f",$hsp->get_aln->percentage_identity), 98.31);
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
is($result->query_length, 2001);
is($result->get_parameter('matrix'), 'blosum62');

is($result->get_statistic('lambda'), 0.318);
is($result->get_statistic('kappa'), 0.135);
is($result->get_statistic('entropy'),0.401 );
is($result->get_statistic('dbentries'), 400);

@valid = ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '6.4e-70', 318],
	   [ 'gi|2367383|gb|AE000509.1|AE000509', 10589, 'AE000509', '0.9992', 59]
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
is($hit->name, 'PIR2:S44629');
is($hit->length, 628);
is($hit->accession, 'PIR2:S44629');
is($hit->significance, '2e-08' );

TODO: {
    local $TODO = 'Raw score parsing broken for GCG-BLAST Hits -- see HSP';
	is($hit->raw_score, 57 );
}

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
my @dcompare = ( ['Contig3700', 5631, 785, '0.0', 785, '0.0', 396, 639, 12, 
		  8723,9434, 1, 4083, 4794, -1],
                 ['Contig3997', 12734, 664, '0.0', 664, '0.0', 335, 401, 0, 
		  1282, 1704, 1, 1546, 1968,-1 ],
                 ['Contig634', 858, 486, '1e-136', 486, '1e-136', 245, 304, 3, 
		  7620, 7941, 1, 1, 321, -1],
                 ['Contig1853', 2314, 339, '1e-91',339, '1e-91', 171, 204, 0,
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
is($r->num_hits, 7);
$hit = $r->next_hit;
is($hit->name, 'gnl|CDD|3919');
is($hit->significance, 0.064);
is($hit->raw_score, 28);
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
my @eqset = qw( 
		
		c200-vs-yeast.BLASTN.m9);
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

# Test Blast parsing with B=0 (WU-BLAST)
$searchio = Bio::SearchIO->new(-file   => test_input_file('no_hsps.blastp'),
			      -format => 'blast');
$result = $searchio->next_result;
is($result->query_name, 'mgri:MG00189.3');
$hit = $result->next_hit;
is($hit->name, 'mgri:MG00189.3');
is($hit->description, 'hypothetical protein 6892 8867 +');
is($hit->score, 3098);
is($hit->significance, '0.');

$hit = $result->next_hit;
is($hit->name, 'fgram:FG01141.1');
is($hit->description, 'hypothetical protein 47007 48803 -');
is($hit->score, 2182);
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


@valid = ( [ 'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', '6e-059', 236],
	      [ 'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', '4e-026', 127],
	      [ 'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', '1e-023', 119]);
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
is($result->query_length, 193);
@hits = $result->hits;
is(scalar(@hits), 10);
is($hits[1]->accession,'1W30');
is($hits[4]->significance,'2e-72');
is($hits[7]->score,'254');
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
is($hits[4]->significance, '0.0');
is($hits[7]->score, 624);

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
is $hit->raw_score, 23;
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
is $hit->raw_score, 120;
is $hit->significance, 3e-27;
$hit = $result->next_hit;
is $hit->name, 'ENSP00000327738';
is $hit->length, 468;
is $hit->accession, 'ENSP00000327738';
is $hit->description, 'pep:known-ccds chromosome:NCBI36:4:189297592:189305643:1 gene:ENSG00000184108 transcript:ENST00000332517 CCDS3851.1';
is $hit->raw_score, 115;
is $hit->significance, 8e-26;

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
