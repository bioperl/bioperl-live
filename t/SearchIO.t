# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;

use vars qw($SKIPXML $LASTXMLTEST); 
use strict;
use lib '.';
use Dumpvalue();
my $dumper = new Dumpvalue();

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use vars qw($NTESTS);
	$NTESTS = 1337;
	$LASTXMLTEST = 67;
	$error = 0;

	use Test;
	plan tests => $NTESTS; 

	eval { require XML::Parser::PerlSAX; 
			 require HTML::Entities; };
	if( $@ ) {
		$SKIPXML = 1;
		print STDERR "XML::Parser::PerlSAX or HTML::Entities not loaded. This means SearchIO::blastxml test cannot be executed. Skipping\n";
		foreach ( 1..$LASTXMLTEST ) {
			skip('No XML::Parser::PerlSAX or HTML::Entities loaded',1);
		}
	}
}

if( $error == 1 ) {
    exit(0);
}


use Bio::SearchIO;
use Bio::Root::IO;
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;

END { 
    unlink 'searchio.out';
    unlink 'searchio.html';
}

ok(1);
my ($searchio, $result,$iter,$hit,$hsp);

if( ! $SKIPXML ) {
    # test with RPSBLAST data first 
    $searchio = new Bio::SearchIO ('-tempfile' => 1,
				   '-format' => 'blastxml',
				   '-file'   => Bio::Root::IO->catfile('t','data','ecoli_domains.rps.xml'));
    
    $result = $searchio->next_result;
    ok($result);    
    ok($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam');
    ok($result->query_name,'gi|1786182|gb|AAC73112.1|');
    ok($result->query_description, '(AE000111) thr operon leader peptide [Escherichia coli]');
    ok($result->query_accession, 'AAC73112.1');
    ok($result->query_length, 21);
    ok($result->algorithm, 'BLASTP');
    ok($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');

    ok($result->available_parameters, 8);
    ok($result->get_parameter('gapext'), 1);
    ok($result->available_statistics, 5);
    ok($result->get_statistic('lambda'), 0.267);

# this result actually has a hit
    $result = $searchio->next_result;
    $hit = $result->next_hit;
    ok($hit->name, 'gnl|Pfam|pfam00742');
    ok($hit->description(), 'HomoS_dh, HomoS dehydrogenase');
    ok($hit->accession, 'pfam00742');
    ok($hit->length, 310);

    $hsp = $hit->next_hsp;
    ok($hsp->query->seq_id, $result->query_name,'query name on HSP');
    ok($hsp->query->seqdesc, $result->query_description,'query desc on HSP');
    ok($hsp->hit->seq_id, $hit->name,'hitname');
    ok($hsp->hit->seqdesc, $hit->description,'hitdesc');
    ok($hsp->pvalue, undef);
    ok(sprintf("%g",$hsp->evalue), sprintf("%g",'1.46134e-90'));
    ok($hsp->score, 838);
    ok($hsp->bits,327.405);
    ok($hsp->query->start, 498);
    ok($hsp->query->end,815);
    ok($hsp->hit->start, 3);
    ok($hsp->hit->end, 310);
    ok($hsp->query->frame,0);
    ok($hsp->hit->frame,0);
    ok(sprintf("%.2f", $hsp->percent_identity), 37.73);
    ok(sprintf("%.4f", $hsp->frac_identical('hit')), 0.3994);
    ok(sprintf("%.4f", $hsp->frac_identical('query')), 0.3868);
    ok(sprintf("%.4f",$hsp->query->frac_identical), 0.3868);

    ok(sprintf("%.4f",$hsp->frac_conserved('total')),0.5245);
    ok(sprintf("%.4f",$hsp->frac_conserved('hit')),0.5552);
    ok(sprintf("%.4f",$hsp->frac_conserved('query')),0.5377);
    ok($hsp->gaps('total'), 26);
    ok($hsp->length('hsp'), 326);
    ok($hsp->query_string, 'LRVCGVANSKALLTNVHGLNLENWQEELAQAKEPF-NLGRLIRLVKEYHLLN----PVIVDCTSSQAVAD-QYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE-GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARET-GRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS');
    ok($hsp->hit_string, 'GVVTGITDSREMLLSRIGLPLEIWKVALRDLEKPRKDLGKLDLTDDAFAVVDDPDIDVVVELTGGIEVARELYLDALEEGKHVVTANKALNASHGDEYLAL---AEKSGVDVLYEAAVAGGIPIIKTLRELLATGDRILKIEGIFNGTTNFILSEMDEKGLPFSDVLAEAQELGYTEADPRDDVEGIDAARKLAILARIAFGIELELDDVYVEGISPITAEDISSADEFGYTLKLLDEAMRQRVEDAESGGEVLRYPTLIPE-------------DHPLASVKGSDNAVAVEGEAYG--PLMFYGPGAGAEPTASAVVADIVRIAR');
    ok($hsp->homology_string, '  V G+ +S+ +L +  GL LE W+  L   ++P  +LG+L      + +++     V+V+ T    VA   Y D L EG HVVT NK  N S  D Y  L   AEKS    LY+  V  G+P+I+ L+ LL  GD ++K  GI +G+ ++I  ++DE G+ FS+    A+E+GYTE DPRDD+ G+D ARKL ILAR   G ELEL D+ +E + P           F   L  LD+    RV  A   G+VLRY   I E             + PL  VK  +NA+A     Y   PL+  G GAG + TA+ V AD++R   ');
    ok(join(' ', $hsp->seq_inds('query', 'gap',1)), '532 548-551 562 649 690');
# one more 
    $hit = $result->next_hit;
    ok($hit);
    
    my $results_left = 8;
    while( $result = $searchio->next_result ) { ok($result); $results_left--; }
    ok $results_left, 0;


    $searchio = new Bio::SearchIO(-format => 'blastxml', 
				  -file => Bio::Root::IO->catfile('t','data','plague_yeast.bls.xml'));

    $result = $searchio->next_result;

    ok($result->database_name, 'yeast.aa');
    ok($result->query_name, 'gi|5763811|emb|CAB53164.1|');
    ok($result->query_description,  'putative transposase [Yersinia pestis]');
    ok($result->query_accession, 'CAB53164.1');
    ok($result->query_length, 340);

    $hit = $result->next_hit;
    ok(! $hit);

    $searchio = new Bio::SearchIO(-format => 'blastxml', 
				  -file => Bio::Root::IO->catfile('t','data','mus.bls.xml'));

    $result = $searchio->next_result;

    ok($result->database_name,'Hs15_up1000');
    ok($result->query_name,'NM_011441_up_1000_chr1_4505586_r');
    ok($result->query_description,'chr1:4505586-4506585');
    ok($result->query_accession,'NM_011441_up_1000_chr1_4505586_r');
    ok($result->query_length,'1000');
    $hit = $result->next_hit;
    ok($hit->name,'NM_001938_up_1000_chr1_93161154_f');
    ok($hit->description,'chr1:93161154-93162153');
    ok($hit->accession,'3153');
    ok($hit->length,'1000');
}
$searchio = new Bio::SearchIO ('-format' => 'blast',
				  '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'));

$result = $searchio->next_result;
# $dumper->dumpValue($result);

ok($result->database_name, 'ecoli.aa');
ok($result->database_entries, 4289);
ok($result->database_letters, 1358990);

ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, qr/^2\.1\.3/);
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_length, 820);
ok($result->get_statistic('kappa'), '0.135');
ok($result->get_statistic('kappa_gapped'), '0.0410');
ok($result->get_statistic('lambda'), '0.319');
ok($result->get_statistic('lambda_gapped'), '0.267');
ok($result->get_statistic('entropy'), '0.383');
ok($result->get_statistic('entropy_gapped'), '0.140');

ok($result->get_statistic('dbletters'), 1358990);
ok($result->get_statistic('dbentries'), 4289);
ok($result->get_statistic('effective_hsplength'), 47);
ok($result->get_statistic('effectivespace'), 894675611);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_parameter('gapopen'), 11);
ok($result->get_parameter('gapext'), 1);
ok($result->get_statistic('S2'), '92');
ok($result->get_statistic('S2_bits'), '40.0');
ok($result->get_parameter('expect'), '1.0e-03');
ok($result->get_statistic('num_extensions'), '82424');


my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 1567],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922', '1e-91', 332],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994', '3e-47', 184]);
my $count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1);
            ok($hsp->query->end, 820);
            ok($hsp->hit->start, 1);
            ok($hsp->hit->end, 820);
            ok($hsp->length('hsp'), 820);
            ok($hsp->start('hit'), $hsp->hit->start);
            ok($hsp->end('query'), $hsp->query->end);
            ok($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            ok($hsp->evalue == 0.0);
            ok($hsp->score, 4058);
            ok($hsp->bits,1567);	    	    
            ok(sprintf("%.2f",$hsp->percent_identity), 98.29);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.9829);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9829);
            ok($hsp->gaps, 0);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.wublastp'));

$result = $searchio->next_result;

ok($result->database_name, 'ecoli.aa');
ok($result->database_letters, 1358990);
ok($result->database_entries, 4289);
ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, qr/^2\.0MP\-WashU/);
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_accession, 'AAC73113.1');

ok($result->query_length, 820);
ok($result->get_statistic('kappa'), 0.136);
ok($result->get_statistic('lambda'), 0.319);
ok($result->get_statistic('entropy'), 0.384);
ok($result->get_statistic('dbletters'), 1358990);
ok($result->get_statistic('dbentries'), 4289);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_statistic('Frame+0_lambda_used'), '0.319');
ok($result->get_statistic('Frame+0_kappa_used'), '0.136');
ok($result->get_statistic('Frame+0_entropy_used'), '0.384');

ok($result->get_statistic('Frame+0_lambda_computed'), '0.319');
ok($result->get_statistic('Frame+0_kappa_computed'), '0.136');
ok($result->get_statistic('Frame+0_entropy_computed'), '0.384');

ok($result->get_statistic('Frame+0_lambda_gapped'), '0.244');
ok($result->get_statistic('Frame+0_kappa_gapped'), '0.0300');
ok($result->get_statistic('Frame+0_entropy_gapped'), '0.180');

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
        ok($qcontigs->[0]->{'start'}, 5);
        ok($qcontigs->[0]->{'stop'}, 812);
        ok($qcontigs->[0]->{'iden'}, 250);
        ok($qcontigs->[0]->{'cons'}, 413);
        # Subject contigs
        ok($scontigs->[0]->{'start'}, 16);
        ok($scontigs->[0]->{'stop'}, 805);
        ok($scontigs->[0]->{'iden'}, 248);
        ok($scontigs->[0]->{'cons'}, 410);
    }

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1);
            ok($hsp->query->end, 820);
            ok($hsp->hit->start, 1);
            ok($hsp->hit->end, 820);
            ok($hsp->length('hsp'), 820);
            
            ok($hsp->evalue == 0.0);
            ok($hsp->pvalue == 0.0);
            ok($hsp->score, 4141);
            ok($hsp->bits,1462.8);	    	    
            ok($hsp->percent_identity, 100);
            ok($hsp->frac_identical('query'), 1.00);
            ok($hsp->frac_identical('hit'), 1.00);
            ok($hsp->gaps, 0);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

# test that add hit really works properly for BLAST objects
# bug 1611
my @hits = $result->hits;
$result->add_hit($hits[0]);
ok($result->num_hits, @hits + 1);

# test WU-BLAST -noseqs option
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.noseqs.wublastp'));

$result = $searchio->next_result;

ok($result->database_name, 'ecoli.aa');
ok($result->database_letters, 1358990);
ok($result->database_entries, 4289);
ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, qr/^2\.0MP\-WashU/);
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_accession, 'AAC73113.1');

ok($result->query_length, 820);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('lambda'), 0.319);
ok($result->get_statistic('entropy'), 0.384);
ok($result->get_statistic('dbletters'), 1358990);
ok($result->get_statistic('dbentries'), 4289);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_statistic('Frame+0_lambda_used'), '0.319');
ok($result->get_statistic('Frame+0_kappa_used'), '0.135');
ok($result->get_statistic('Frame+0_entropy_used'), '0.384');

ok($result->get_statistic('Frame+0_lambda_computed'), '0.319');
ok($result->get_statistic('Frame+0_kappa_computed'), '0.135');
ok($result->get_statistic('Frame+0_entropy_computed'), '0.384');

ok($result->get_statistic('Frame+0_lambda_gapped'), '0.244');
ok($result->get_statistic('Frame+0_kappa_gapped'), '0.0300');
ok($result->get_statistic('Frame+0_entropy_gapped'), '0.180');

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113', '0', 4141],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922', '6.6e-93', 907],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994', '2.8e-47', 483]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1);
            ok($hsp->query->end, 820);
            ok($hsp->hit->start, 1);
            ok($hsp->hit->end, 820);
            ok($hsp->length('hsp'), 820);
            
            ok($hsp->evalue == 0.0);
            ok($hsp->pvalue == 0.0);
            ok($hsp->score, 4141);
            ok($hsp->bits,1462.8);	    	    
            ok($hsp->percent_identity, 100);
            ok($hsp->frac_identical('query'), 1.00);
            ok($hsp->frac_identical('hit'), 1.00);
            ok($hsp->gaps, 0);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

# test tblastx 
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

$result = $searchio->next_result;
ok($result->database_name, 'ecoli.nt');
ok($result->database_letters, 4662239);
ok($result->database_entries, 400);
ok($result->algorithm, 'TBLASTX');
ok($result->algorithm_version, qr/^2\.1\.2/);
ok($result->query_name, 'HUMBETGLOA');
ok($result->query_description, 'Human haplotype C4 beta-globin gene, complete cds.');
ok($result->query_length, 3002);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('lambda'), 0.318);
ok($result->get_statistic('entropy'), 0.401);
ok($result->get_statistic('dbletters'), 4662239);
ok($result->get_statistic('dbentries'), 400);
ok($result->get_statistic('T'), 13);
ok($result->get_statistic('X1'), 16);
ok($result->get_statistic('X1_bits'), 7.3);
ok($result->get_statistic('X2'), 0);
ok($result->get_statistic('X2_bits'), '0.0');
ok($result->get_statistic('S1'), 41);
ok($result->get_statistic('S1_bits'), 21.7);
ok($result->get_statistic('S2'), 53);
ok($result->get_statistic('S2_bits'), 27.2);

ok($result->get_statistic('decayconst'), 0.1);

ok($result->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AE000479.1|AE000479', 10934, 'AE000479', '0.13', 34],
	   [ 'gb|AE000302.1|AE000302', 10264, 'AE000302', '0.61', 31],
	   [ 'gb|AE000277.1|AE000277', 11653, 'AE000277', '0.84', 31]);
$count = 0;

while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1057);
            ok($hsp->query->end, 1134);
            ok($hsp->query->strand, 1);
            ok($hsp->strand('query'), $hsp->query->strand);
            ok($hsp->hit->end, 5893);
            ok($hsp->hit->start, 5816);
            ok($hsp->hit->strand, -1);
            ok($hsp->strand('sbjct'), $hsp->subject->strand);
            ok($hsp->length('hsp'), 26);
            
            ok($hsp->evalue == 0.13);
            ok($hsp->score, 67);
            ok($hsp->bits,33.6);
            ok(sprintf("%.2f",$hsp->percent_identity), 42.31);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), '0.4231');
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.4231');
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 1);
            ok($hsp->gaps, 0);	    
            ok($hsp->query_string, 'SAYWSIFPPLGCWWSTLGPRGSLSPL');
            ok($hsp->hit_string, 'AAVWALFPPVGSQWGCLASQWRTSPL');
            ok($hsp->homology_string, '+A W++FPP+G  W  L  +   SPL');
            ok(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '355 364 365 367 368 370 371 373-375');
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => File::Spec->catfile(qw(t data HUMBETGLOA.FASTA)) );
$result = $searchio->next_result;
ok($result->database_name, qr/dros_clones.2.5/);
ok($result->database_letters, 112936249);
ok($result->database_entries, 657);
ok($result->algorithm, 'FASTN');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, "HUMBETGLOA");
ok($result->query_description, "Human haplotype C4 beta-globin gene, complete cds.");
ok($result->query_length, 3002);
ok($result->get_parameter('gapopen'), -16);
ok($result->get_parameter('gapext'), -4);
ok($result->get_parameter('ktup'), 6);

ok($result->get_statistic('lambda'), 0.0823);
ok($result->get_statistic('dbletters'), 112936249);
ok($result->get_statistic('dbentries'), 657);

@valid = ( [ 'BACR21I23', 73982, 'BACR21I23', '0.017', 44.2],
	   [ 'BACR40P19', 73982, 'BACR40P19', '0.017', 44.2],
	   [ 'BACR30L17', 32481, 'BACR30L17', '0.018', 44.1]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );
    ok($hit->rank, $count + 1);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 31);
            ok($hsp->query->end, 289);
            ok($hsp->query->strand, -1);
            ok($hsp->hit->end, 65167);
            ok($hsp->hit->start, 64902);
            ok($hsp->hit->strand, 1);
            ok($hsp->length('hsp'), 267);	    
            ok($hsp->evalue == 0.017);
            ok($hsp->score, 134.5);
            ok($hsp->bits,44.2);
            ok(sprintf("%.2f",$hsp->percent_identity), '57.30');
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.5907); 
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5752); 
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps, 159);
            ok($hsp->gaps('query'), 8);
            ok($hsp->gaps('hit'),1);
            ok($hsp->query_string, 'GATTAAAACCTTCTGGTAAGAAAAGAAAAAATATATATATATATATATGTGTATATGTACACACATACATATACATATATATGCATTCATTTGTTGTTGTTTTTCTTAATTTGCTCATGCATGCTA----ATAAATTATGTCTAAAAATAGAAT---AAATACAAATCAATGTGCTCTGTGCATTA-GTTACTTATTAGGTTTTGGGAAACAAGAGGTAAAAAACTAGAGACCTCTTAATGCAGTCAAAAATACAAATAAATAAAAAGTCACTTACAACCCAAAGTGTGACTATCAATGGGGTAATCAGTGGTGTCAAATAGGAGGT');
            ok($hsp->hit_string, 'GATGTCCTTGGTGGATTATGGTGTTAGGGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATATAAAATATAATATAAAATATAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATAT-AATATAAAATATAAAATAAAATATAATATAAAATATAATATAAAATATAATATAAAATATAATATAAAATA');
            ok($hsp->homology_string, '                              :::::::::::::::::: : ::::: :: : : ::: ::::: ::::::::  ::  :: : :   : : : : :  ::    : :: ::   ::    : ::: :::     :::::: :::   ::::: ::  :::  :    :    : ::   :::  : ::   : :   : : :: :   :: : : :: : :       ::  : : ::: ::: ::  ::::: ::: : :  :: ::   ::: : : : ::: ::   '.' 'x60);
            ok(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '33 37 39 41 43 47-49 52 55 56 58 60 64 70 71 74 78 82 84 86 87 90-96 98 100 103 105 107 110-112 114 117 119 121-123 125 127-129 132 134 135 139-141 143 145-148 150-153 155 156 160 161 164 170 173 180-184 188 192 194 196-198 201 204 206-209 212 213 215 217 219 221 223-225 227 229 232 233 236 237 246 252 256 258 260 263 269 271');
            ok(join(' ', $hsp->seq_inds('query', 'conserved',1)), '31 32 34-36 38 40 42 44-46 50 51 53 54 57 59 61-63 65-69 72 73 75-77 79-81 83 85 88 89 97 99 101 102 104 106 108 109 113 115 116 118 120 124 126 130 131 133 136-138 141 142 144 149 154 157-159 162 163 165-172 174-179 185-187 189-191 193-195 199 200 202 203 205 210 211 214 216 218 220 222 226 228 230 231 234 235 238-245 247-251 253-255 257 259 261 262 264-268 270 272-289');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            ok(sprintf("%.2f",$hsp->get_aln->percentage_identity()), '59.30');
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => File::Spec->catfile(qw(t data cysprot1.FASTA)));
$result = $searchio->next_result;
ok($result->database_name, qr/ecoli.aa/);
ok($result->database_letters, 1358987);
ok($result->database_entries, 4289);
ok($result->algorithm, 'FASTP');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, 'CYS1_DICDI');
ok($result->query_length, 343);
ok($result->get_parameter('gapopen'), -12);
ok($result->get_parameter('gapext'), -2);
ok($result->get_parameter('ktup'), 2);

ok($result->get_statistic('lambda'), 0.1456);
ok($result->get_statistic('dbletters'), 1358987);
ok($result->get_statistic('dbentries'), 4289);


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309', 1.2, 29.2],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148', 2.1, 27.4],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494',  2.1, 25.9]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 125);
            ok($hsp->query->end, 305);
            ok($hsp->query->strand, 0);
            ok($hsp->hit->start, 2);
            ok($hsp->hit->end, 181);
            ok($hsp->hit->strand, 0);
            ok($hsp->length('hsp'), 188);	    
            ok($hsp->evalue == 1.2);
            ok($hsp->score, 109.2);
            ok($hsp->bits,29.2);
            ok(sprintf("%.2f",$hsp->percent_identity), 23.94);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.2486);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.2500');
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps('query'), 7);
            ok($hsp->gaps, 49);	    
            ok($hsp->query_string, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
            ok($hsp->hit_string, (' 'x29).'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
            ok($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                      '.' 'x60);
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            ok(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 26.01);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

ok($result->hits, 8);
$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => File::Spec->catfile(qw(t data cysprot_vs_gadfly.FASTA)) );
$result = $searchio->next_result;
ok($result->database_name, qr/gadflypep2/);
ok($result->database_letters, 7177762);
ok($result->database_entries, 14334);
ok($result->algorithm, 'FASTP');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, 'cysprot.fa');
ok($result->query_length, 2385);
ok($result->get_parameter('gapopen'), -12);
ok($result->get_parameter('gapext'), -2);
ok($result->get_parameter('ktup'), 2);
ok($result->get_parameter('matrix'), 'BL50');

ok($result->get_statistic('lambda'), 0.1397);
ok($result->get_statistic('dbletters'), 7177762 );
ok($result->get_statistic('dbentries'), 14334);


@valid = ( [ 'Cp1|FBgn0013770|pp-CT20780|FBan0006692', 341, 
	     'FBan0006692', '3.1e-59', 227.8],
	   [ 'CG11459|FBgn0037396|pp-CT28891|FBan0011459', 336, 
	     'FBan0011459', '6.4e-41',  166.9],
	   [ 'CG4847|FBgn0034229|pp-CT15577|FBan0004847', 390, 
	     'FBan0004847',  '2.5e-40', 165.2]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1373);
            ok($hsp->query->end, 1706);
            ok($hsp->query->strand, 0);
            ok($hsp->hit->start, 5);
            ok($hsp->hit->end, 341);
            ok($hsp->hit->strand, 0);
            ok($hsp->length('hsp'), 345);	    
            ok(sprintf("%g",$hsp->evalue), sprintf("%g",'3.1e-59') );
            ok($hsp->score, 1170.6);
            ok($hsp->bits,227.8);
            ok(sprintf("%.2f",$hsp->percent_identity), 53.04);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.5479);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.5430');
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps('query'), 11);
            ok($hsp->gaps, 194);
            ok($hsp->hit_string, (' 'x26).'MRTAVLLPLLAL----LAVAQA-VSFADVVMEEWHTFKLEHRKNYQDETEERFRLKIFNENKHKIAKHNQRFAEGKVSFKLAVNKYADLLHHEFRQLMNGFNYTLHKQLRAADESFKGVTFISPAHVTLPKSVDWRTKGAVTAVKDQGHCGSCWAFSSTGALEGQHFRKSGVLVSLSEQNLVDCSTKYGNNGCNGGLMDNAFRYIKDNGGIDTEKSYPYEAIDDSCHFNKGTVGATDRGFTDIPQGDEKKMAEAVATVGPVSVAIDASHESFQFYSEGVYNEPQCDAQNLDHGVLVVGFGTDESGED---YWLVKNSWGTTWGDKGFIKMLRNKENQCGIASASSYPLV');
            ok($hsp->query_string, 'SNWGNNGYFLIERGKNMCGLAACASYPIPQVMNPTLILAAFCLGIASATLTFDHSLEAQWTKWKAMHNRLY-GMNEEGWRRAVWEKNMKMIELHNQEYREGKHSFTMAMNAFGDMTSEEFRQVMNGFQ---NRKPR------KGKVFQEPLFYEAPRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQDNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPK-QEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCSSEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTVMTPLLLLAVLCLGTALATPKFDQTFNAQWHQWKSTHRRLYGTNEE');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            ok(sprintf("%.2f",$hsp->get_aln->percentage_identity()), 56.13);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;
ok($result->hits, 21);

# test on TFASTXY
$searchio = new Bio::SearchIO(-format => 'fasta',
			      -file   => File::Spec->catfile(qw(t data 5X_1895.FASTXY)));
$result = $searchio->next_result;
ok($result->database_name, qr/yeast_nrpep.fasta/);
ok($result->database_letters, 4215311);
ok($result->database_entries, 9190);
ok($result->algorithm, 'FASTY');
ok($result->algorithm_version, '3.4t07');
ok($result->query_name, '5X_1895.fa');
ok($result->query_length, 7972);
ok($result->get_parameter('gapopen'), -14);
ok($result->get_parameter('gapext'), -2);
ok($result->get_parameter('ktup'), 2);
ok($result->get_parameter('matrix'), 'BL50');

ok($result->get_statistic('lambda'), 0.1711);
ok($result->get_statistic('dbletters'), 4215311);
ok($result->get_statistic('dbentries'), 9190);


@valid = ( [ 'NR_SC:SW-YNN2_YEAST', 1056, 'NR_SC:SW-YNN2_YEAST','1.6e-154', '547.0'],
	   [ 'NR_SC:SW-MPCP_YEAST', 311, 'NR_SC:SW-MPCP_YEAST', '1.3e-25', 117.1],
	   [ 'NR_SC:SW-YEO3_YEAST', 300, 'NR_SC:SW-YEO3_YEAST', '5.7e-05', 48.5]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 2180);
            ok($hsp->query->end, 5623);
            ok($hsp->query->strand, 1);
            ok($hsp->hit->start, 3);
            ok($hsp->hit->end, 1053);
            ok($hsp->hit->strand, 0);
            ok($hsp->length('hsp'), 1165);
            
            ok(sprintf("%g",$hsp->evalue), sprintf("%g",'1.6e-154'));
            ok($hsp->score, 2877.6);
            ok($hsp->bits,'547.0');
            ok(sprintf("%.2f",$hsp->percent_identity), 51.67);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.5244);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5728);
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps, 678);	    
            ok($hsp->query_string, 'RKQLDPRIPALINNGVKANHRSFFVMVGDKGRDQVCPGMQAAMRFD*HRCR/LVNLHFLLSQARVSSRPSVLWCYKKD-LGFTT*VAASENLQQTIYFRPIATSHRKKREAKIKRDVKRGIRDANEQDPFELFVTVTDIRYTYYKDSAKILGQTFGMLVLQDYEAITPNLLARTIETVEGGGIVVLLLKTMSSLKQLYAMAM/DKL*CRDGVE*SDFS*LLI*DVHSRYRTDAHQFVQPRFNERFILSLGSNPDCLVLDDELNVLPLSKGKDIQIGKAGEEDDRGRKRKAEELKEMKENLEGVDIVGSLAKLAKTVDQAKAILTFVEAISEKNLSSTVALTAGRGRGKSAALGLAIGAALAHDYSNIFVTSPDPENLKTLFEFVFKALDALGYEEHIDYDVVQSTNPDFKKAIVRVNIFRGHRQTIQYISPEDSHVLGQAELVIIDEAAAIPLPLVRKLIGPYLVFMASTINGYEGTGRSLSIKLIQQLREQTRPSITKDSENAAASSAGSSSKAAAAGRSGAGLVRSLREIKLDEPIRYSPGDNVEKWLNNLLCLDATIVSK---SIQGCPHPSKCELYYVNRDTLFSYHPASEVFLQRMMALYVASHYKNSPNDLQMLSDAPAHHLFVLLPPIDEND-NTLPDPLVVLQVALEGNISREAILKEMAQSGMRSSGDMIPWIISTQFQDNDFATLSGARVVRIATHPDYARMGYGSRAMEALESFYNGTSYNFDDVPVDMGESFAD\VPRSDL*VTSFIPFPQNRTSTECVSQNANLQNDTIAIRDPSRMPPLLQRLSERKPETLDYLGVSFGLTRDLLRFWKKGGFTPLYASQKENALTGEYTFVMLKVLASAGGGGEWLGAFAQGMSCLLLQDEVHMGND*RL*TDFRQRFMNLLSYEAFKKFDASIALSILESTVPRNSPSPAP----KLLTNTELSSLLTPFDIKRLESYADSMLDYHVVLDLVPTIASLFFGKRLETS--LPPAQQAILLALGLQRKNVEALENELGITSTQTLALFGKVLRKMTKSLEDIRKASIASELP-----AEPTLAGRSANGSNKFVALQQTIEQDLADSAVQLNGEDDDASKKEQRELLNTLNMEEFAI-DQGGDWTEAEKQVERLASGKGGTRLSSTVSVKVDKLDD\AKRRRRRARMRVPRMRRR');
            ok($hsp->hit_string, 'KKAIDSRIPSLIRNGVQTKQRSIFVIVGDRARNQ------------------LPNLHYLMMSADLKMNKSVLWAYKKKLLGFT--------------------SHRKKRENKIKKEIKRGTREVNEMDPFESFISNQNIRYVYYKESEKILGNTYGMCILQDFEALTPNLLARTIETVEGGGIVVILLKSMSSLKQLYTMTM-D--------------------VHARYRTEAHGDVVARFNERFILSLGSNPNCLVVDDELNVLPLSGAKNVKPLPPKEDDELPPKQL--ELQELKESLEDVQPAGSLVSLSKTVNQAHAILSFIDAISEKTLNFTVALTAGRGRGKSAALGISIAAAVSHGYSNIFVTSPSPENLKTLFEFIFKGFDALGYQEHIDYDIIQSTNPDFNKAIVRVDIKRDHRQTIQYIVPQDHQVLGQAELVVIDEAAAIPLPIVKNLLGPYLVFMASTINGYEGTGRSLSLKLIQQLRNQNNTSGRESTQTAVVSRDNKEKDSHLHSQS-----RQLREISLDEPIRYAPGDPIEKWLNKLLCLDVTLIKNPRFATRGTPHPSQCNLFVVNRDTLFSYHPVSENFLEKMMALYVSSHYKNSPNDLQLMSDAPAHKLFVLLPPIDPKDGGRIPDPLCVIQIALEGEISKESVRNSLSR-GQRAGGDLIPWLISQQFQDEEFASLSGARIVRIATNPEYASMGYGSRAIELLRDYFEGKF-------TDMSE---D-VRPKDYSI--------KRVSDKELAKT-NLLKDDVKLRDAKTLPPLLLKLSEQPPHYLHYLGVSYGLTQSLHKFWKNNSFVPVYLRQTANDLTGEHTCVMLNVLE--GRESNWLVEFAK---------------------DFRKRFLSLLSYD-FHKFTAVQALSVIESSKKAQDLSDDEKHDNKELTRTHLDDIFSPFDLKRLDSYSNNLLDYHVIGDMIPMLALLYFGDKMGDSVKLSSVQSAILLAIGLQRKNIDTIAKELNLPSNQTIAMFAKIMRKMSQYFRQLLSQSIEETLPNIKDDAIAEMDGEEIKNYNAAEALDQ-MEEDLEEAG----SEAVQAMREKQKELINSLNLDKYAINDNSEEWAESQKSLEIAAKAKGVVSLKTGKKRTTEKAED-IYRQEMKA-MKKPRKSKK');
            ok($hsp->homology_string, '.: .: :::.:: :::....::.::.:::..:.:                  : :::.:. .: ..   :::: :::  ::::                    ::::::: :::...::: :..::.:::: :..  .:::.:::.: ::::.:.:: .:::.::.:::::::::::::::::::.:::.::::::::.:.: :                    ::.::::.::  :  ::::::::::::::.:::.:::::::::: .:...     :.:.   :.   ::.:.::.:: :. .:::..:.:::.::.:::.:..:::::.:. :::::::::::::::::..:.::..: :::::::::.::::::::::.::..:::::.::::::..:::::::.::::::.: : :::::::: :.: .::::::::.::::::::::.:..:.::::::::::::::::::::::.:::::::.:.  :  .....:..:  .. . .   ..:     :.::::.:::::::.::: .:::::.:::::.:....   . .: ::::.:.:. :::::::::::.:: ::..::::::.:::::::::::..::::::.::::::::: .: . .:::: :.:.::::.::.:.. . ... :.:..::.:::.:: ::::..::.:::::.:::::.:.:: :::::::.: :.....:         .::.:   : :  .:  .        .:.: . .... :: .: . .:: . .:::: .:::. :. : :::::.:::..: .:::...:.:.:  :  : ::::.: :::.::   :  ..::  ::.                     :::.::..::::. :.:: :  :::..::.   .. :       : :: :.:.....:::.:::.::....:::::. :..: .: :.:: ..  :  :  .:.:::::.::::::.... .::.. :.::.:.:.:..:::.. .... . ::   ::     :   . :.  .. :   ::.: .:.:: ...    .:  .: ...:.::.:.::....:: :.. .:.:..:..:  :..:: . :..  .  ..: .:   :.. .: :. ::  ..');
            # note: the reason this is not the same percent id above
            # is we are calculating average percent id
            ok(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity()),
               '51.77');
            ok(sprintf("%.2f",$hsp->get_aln->average_percentage_identity()),
               '58.41');
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;
ok($result->hits, 58);
# test for MarkW bug in blastN

$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','a_thaliana.blastn'));


$result = $searchio->next_result;
ok($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS, GSS,or phase 0, 1 or 2 HTGS sequences) ');
ok($result->database_letters, 4677375331);
ok($result->database_entries, 1083200);
ok($result->algorithm, 'BLASTN');
ok($result->algorithm_version, qr/^2\.2\.1/);
ok($result->query_name, '');
ok($result->query_length, 60);
ok($result->get_parameter('gapopen'), 5);
ok($result->get_parameter('gapext'), 2);
ok($result->get_parameter('ktup'), undef);

ok($result->get_statistic('lambda'), 1.37);
ok($result->get_statistic('kappa'), 0.711);
ok($result->get_statistic('entropy'),1.31 );
ok($result->get_statistic('T'), 0);
ok($result->get_statistic('A'), 30);
ok($result->get_statistic('X1'), '6');
ok($result->get_statistic('X1_bits'), 11.9);
ok($result->get_statistic('X2'), 15);
ok($result->get_statistic('X2_bits'), 29.7);
ok($result->get_statistic('S1'), 12);
ok($result->get_statistic('S1_bits'), 24.3);
ok($result->get_statistic('S2'), 17);
ok($result->get_statistic('S2_bits'), 34.2);

ok($result->get_statistic('dbentries'), 1083200);

@valid = ( [ 'gb|AY052359.1|', 2826, 'AY052359', '3e-18', 96, 1, 60, 
	     '1.0000'],
	   [ 'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 96, 1, 60, 
	     '1.0000' ],
	   [ 'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.04', 42, 35, 55, 
	     '0.3500']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );
    ok($hit->start, shift @$d);
    ok($hit->end,shift @$d);    
    ok(sprintf("%.4f",$hit->frac_aligned_query), shift @$d);
    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1);
            ok($hsp->query->end, 60);
            ok($hsp->query->strand, 1);
            ok($hsp->hit->start, 154);
            ok($hsp->hit->end, 212);
            ok($hsp->hit->strand, 1);
            ok($hsp->length('hsp'), 60);	    
            ok(sprintf("%g",$hsp->evalue), sprintf("%g",'3e-18'));
            ok($hsp->score, 48);
            ok($hsp->bits,95.6);
            ok(sprintf("%.2f",$hsp->percent_identity), 96.67);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.9667);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9831);
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps('query'), 0);
            ok($hsp->gaps('hit'), 1);
            ok($hsp->gaps, 1);	    
            ok($hsp->query_string, 'aggaatgctgtttaattggaatcgtacaatggagaatttgacggaaatagaatcaacgat');
            ok($hsp->hit_string, 'aggaatgctgtttaattggaatca-acaatggagaatttgacggaaatagaatcaacgat');
            ok($hsp->homology_string, '|||||||||||||||||||||||  |||||||||||||||||||||||||||||||||||');
            ok(sprintf("%.2f",$hsp->get_aln->overall_percentage_identity), 96.67);
            ok(sprintf("%.2f",$hsp->get_aln->percentage_identity), 98.31);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
} 
ok @valid, 0;

#WU-BlastX test

$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','dnaEbsub_ecoli.wublastx'));

$result = $searchio->next_result;
ok($result->database_name, 'ecoli.aa');
ok($result->database_letters, 1358990);
ok($result->database_entries, 4289);
ok($result->algorithm, 'BLASTX');
ok($result->algorithm_version, qr/^2\.0MP\-WashU/);
ok($result->query_name, 'gi|142864|gb|M10040.1|BACDNAE');
ok($result->query_description, 'B.subtilis dnaE gene encoding DNA primase, complete cds');
ok($result->query_accession, 'M10040.1');
ok($result->query_length, 2001);
ok($result->get_parameter('matrix'), 'blosum62');

ok($result->get_statistic('lambda'), 0.318);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('entropy'),0.401 );

ok($result->get_statistic('dbentries'), 4289);

@valid = ( [ 'gi|1789447|gb|AAC76102.1|', 581, 'AAC76102', '1.1e-74', 671]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );
    ok(sprintf("%.4f",$hit->frac_identical('query')), '0.3640');
    ok(sprintf("%.4f",$hit->frac_identical('hit')), '0.3660');
    ok(sprintf("%.4f",$hit->frac_conserved('query')), '0.5370');
    ok(sprintf("%.4f",$hit->frac_conserved('hit')), '0.5400');
    ok(sprintf("%.4f",$hit->frac_aligned_query), '0.6200');
    ok(sprintf("%.4f",$hit->frac_aligned_hit), '0.7100');

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 21);
            ok($hsp->query->end, 1265);
            ok($hsp->query->strand, 1);
            ok($hsp->hit->start, 1);
            ok($hsp->hit->end, 413);
            ok($hsp->hit->strand, 0);
            ok($hsp->length('hsp'), 421);	    
            ok(sprintf("%g",$hsp->evalue), sprintf("%g",'1.1e-74'));
            ok(sprintf("%g",$hsp->pvalue), sprintf("%g",'1.1e-74'));
            ok($hsp->score,671);
            ok($hsp->bits,265.8);
            ok(sprintf("%.2f",$hsp->percent_identity), 35.87);
            
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.3639);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.3656);
            ok(sprintf("%.4f",$hsp->frac_conserved('query')), 0.5373);
            ok(sprintf("%.2f",$hsp->frac_conserved('hit')), 0.54);
            
            ok(sprintf("%.4f",$hsp->frac_identical('hsp')), 0.3587);
            ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.5297);
            
            ok($hsp->query->frame(), 2);
            ok($hsp->hit->frame(), 0);
            ok($hsp->gaps('query'), 6);
            ok($hsp->gaps('hit'), 8);
            ok($hsp->gaps, 14);	    
            ok($hsp->query_string, 'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ');
            ok($hsp->hit_string, 'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ');
            ok($hsp->homology_string, 'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q');
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

#Trickier WU-Blast
$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','tricky.wublast'));
$result = $searchio->next_result;
my $hits_left = 1;
while (my $hit = $result->next_hit) {
	# frac_aligned_hit used to be over 1, frac_identical & frac_conserved are still too wrong
	#ok(sprintf("%.3f",$hit->frac_identical), '0.852');
    #ok(sprintf("%.3f",$hit->frac_conserved), '1.599');
    ok(sprintf("%.2f",$hit->frac_aligned_query), '0.92');
    ok(sprintf("%.2f",$hit->frac_aligned_hit), '0.91');
    $hits_left--;
}
ok $hits_left, 0;

# More frac_ method testing, this time on ncbi blastn
$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','frac_problems.blast'));
my @expected = ("1.000", "0.943");
while (my $result = $searchio->next_result) {
    my $hit = $result->next_hit;
    ok $hit->frac_identical, shift @expected;
}
ok @expected, 0;

#WU-TBlastN test

$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','dnaEbsub_ecoli.wutblastn'));

$result = $searchio->next_result;
ok($result->database_name, 'ecoli.nt');
ok($result->database_letters, 4662239);
ok($result->database_entries, 400);
ok($result->algorithm, 'TBLASTN');
ok($result->algorithm_version, qr/^2\.0MP\-WashU/);
ok($result->query_name, 'gi|142865|gb|AAA22406.1|');
ok($result->query_description, 'DNA primase');
ok($result->query_accession, 'AAA22406.1');
ok($result->query_length, 603);
ok($result->get_parameter('matrix'), 'blosum62');

ok($result->get_statistic('lambda'), '0.320');
ok($result->get_statistic('kappa'), 0.136);
ok($result->get_statistic('entropy'),0.387 );

ok($result->get_statistic('dbentries'), 400);

@valid = ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '1.4e-73', 671]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 1);
            ok($hsp->query->end, 415);
            ok($hsp->query->strand, 0);
            ok($hsp->hit->start, 4778);
            ok($hsp->hit->end, 6016);
            ok($hsp->hit->strand, 1);
            ok($hsp->length('hsp'), 421);	    
            ok(sprintf("%g",$hsp->evalue), sprintf("%g",'1.4e-73'));
            ok(sprintf("%g",$hsp->pvalue), sprintf("%g",'1.4e-73'));
            ok($hsp->score,671);
            ok($hsp->bits,265.8);
            ok(sprintf("%.2f",$hsp->percent_identity), 35.87);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.3656);	    
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.3639);
            ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.5297);
            ok($hsp->query->frame(), 0);
            ok($hsp->hit->frame(), 1);
            ok($hsp->gaps('query'), 6);
            ok($hsp->gaps('hit'), 8);
            ok($hsp->gaps, 14);	    
            ok($hsp->query_string, 'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ');
            ok($hsp->hit_string, 'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ');
            ok($hsp->homology_string, 'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q');
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok $count, 1;

# WU-BLAST TBLASTX
$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data','dnaEbsub_ecoli.wutblastx'));

$result = $searchio->next_result;
ok($result->database_name, 'ecoli.nt');
ok($result->database_letters, 4662239);
ok($result->database_entries, 400);
ok($result->algorithm, 'TBLASTX');
ok($result->algorithm_version, qr/^2\.0MP\-WashU/);
ok($result->query_name, 'gi|142864|gb|M10040.1|BACDNAE');
ok($result->query_description, 'B.subtilis dnaE gene encoding DNA primase, complete cds');
ok($result->query_accession, 'M10040.1');
ok($result->query_length, 2001);
ok($result->get_parameter('matrix'), 'blosum62');

ok($result->get_statistic('lambda'), 0.318);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('entropy'),0.401 );
ok($result->get_statistic('dbentries'), 400);

@valid = ( [ 'gi|1789441|gb|AE000388.1|AE000388', 10334, 'AE000388', '6.4e-70', 318],
	   [ 'gi|2367383|gb|AE000509.1|AE000509', 10589, 'AE000509', '0.9992', 59]
	   );
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    # using e here to deal with 0.9992 coming out right here as well
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hspcounter = 0;
        while( my $hsp = $hit->next_hsp ) {
            $hspcounter++;
            if( $hspcounter == 3 ) {
                # let's actually look at the 3rd HSP
                ok($hsp->query->start, 441);
                ok($hsp->query->end, 617);
                ok($hsp->query->strand, 1);
                ok($hsp->hit->start, 5192);
                ok($hsp->hit->end, 5368);
                ok($hsp->hit->strand, 1);
                ok($hsp->length('hsp'), 59);	    
                ok(sprintf("%g",$hsp->evalue), sprintf("%g",'6.4e-70'));
                ok(sprintf("%g",$hsp->pvalue), sprintf("%g",'6.4e-70'));
                ok($hsp->score,85);
                ok($hsp->bits,41.8);
                ok(sprintf("%.2f",$hsp->percent_identity), '32.20');
                ok(sprintf("%.3f",$hsp->frac_identical('hit')), 0.322);
                ok(sprintf("%.3f",$hsp->frac_identical('query')), 0.322);
                ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.4746);
                ok($hsp->query->frame(), 2);
                ok($hsp->hit->frame(), 1);
                ok($hsp->gaps('query'), 0);
                ok($hsp->gaps('hit'), 0);
                ok($hsp->gaps, 0);	    
                ok($hsp->query_string, 'ALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGY');
                ok($hsp->hit_string, 'ARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY');
                ok($hsp->homology_string, 'A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y');
                last;
            }
        }
        ok $hspcounter, 3;
    }
    elsif( $count == 1 ) {
        my $hsps_to_do = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 587);
            ok($hsp->query->end, 706);
            ok($hsp->query->strand, -1);
            ok($hsp->hit->start, 4108);
            ok($hsp->hit->end, 4227);
            ok($hsp->hit->strand, -1);
            ok($hsp->length('hsp'), 40);	    
            ok($hsp->evalue == '7.1');
            ok($hsp->pvalue == '1.00');
            ok($hsp->score,59);
            ok($hsp->bits,29.9);
            ok(sprintf("%.2f",$hsp->percent_identity), '37.50');
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.3750');
            ok(sprintf("%.4f",$hsp->frac_identical('query')), '0.3750');
            ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), '0.4750');
            ok($hsp->query->frame(), 2);
            ok($hsp->hit->frame(), 2);
            ok($hsp->gaps('query'), 0);
            ok($hsp->gaps('hit'), 0);
            ok($hsp->gaps, 0);
            ok($hsp->query_string, 'WLPRALPEKATTAP**SWIGNMTRFLKRSKYPLPSSRLIR');
            ok($hsp->hit_string, 'WLSRTTVGSSTVSPRTFWITRMKVKLSSSKVTLPSTKSTR');
            ok($hsp->homology_string, 'WL R     +T +P   WI  M   L  SK  LPS++  R');
            $hsps_to_do--;
            last;
        }
        ok $hsps_to_do, 0;
    }
    last if( $count++ > @valid );
}
ok $count, 2;

# Do a multiblast report test
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','multi_blast.bls'));

@expected = qw(CATH_RAT CATL_HUMAN CATL_RAT PAPA_CARPA);
my $results_left = 4;
while( my $result = $searchio->next_result ) {
    ok($result->query_name, shift @expected, "Multiblast query test");
    $results_left--;
}
ok $results_left, 0;

# Test GCGBlast parsing

$searchio = new Bio::SearchIO('-format' => 'blast',
			      '-file'   => Bio::Root::IO->catfile('t','data', 'test.gcgblast'));
$result = $searchio->next_result();

ok($result->query_name, '/v0/people/staji002/test.gcg');
ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, '2.2.1 [Apr-13-2001]');
ok($result->database_name, 'pir');
ok($result->database_entries, 274514);
ok($result->database_letters, 93460074);

$hit = $result->next_hit;
ok($hit->name, 'PIR2:S44629');
ok($hit->length, 628);
ok($hit->accession, 'PIR2:S44629');
skip('Significance parsing broken for GCG-BLAST Hits -- see HSP',$hit->significance, '2e-08' );
skip('Raw score parsing broken for GCG-BLAST Hits -- see HSP',$hit->raw_score, 57 );

$hsp = $hit->next_hsp;
ok(sprintf("%g",$hsp->evalue), sprintf("%g",'2e-08'));
ok($hsp->bits, '57.0');
ok($hsp->score, 136);
ok(int($hsp->percent_identity), 28);
ok(sprintf("%.2f",$hsp->frac_identical('query')), 0.29);
ok($hsp->frac_conserved('total'), 69/135);
ok($hsp->gaps('total'), 8);
ok($hsp->gaps('hit'), 6);
ok($hsp->gaps('query'), 2);

ok($hsp->hit->start, 342);
ok($hsp->hit->end, 470);
ok($hsp->query->start, 3);
ok($hsp->query->end, 135);

ok($hsp->query_string, 'CAAEFDFMEKETPLRYTKTXXXXXXXXXXXXXXRKIISDMWGVLAKQQTHVRKHQFDHGELVYHALQLLAYTALGILIMRLKLFLTPYMCVMASLICSRQLFGW--LFCKVHPGAIVFVILAAMSIQGSANLQTQ');
ok($hsp->hit_string, 'CSAEFDFIQYSTIEKLCGTLLIPLALISLVTFVFNFVKNT-NLLWRNSEEIG----ENGEILYNVVQLCCSTVMAFLIMRLKLFMTPHLCIVAALFANSKLLGGDRISKTIRVSALVGVI-AILFYRGIPNIRQQ');
ok($hsp->homology_string, 'C+AEFDF++  T  +   T                 + +   +L +    +     ++GE++Y+ +QL   T +  LIMRLKLF+TP++C++A+L  + +L G   +   +   A+V VI A +  +G  N++ Q');


$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

$result = $searchio->next_result;

ok($result);
$hit = $result->next_hit;
ok($hit->accession, 'AE000479');
ok($hit->bits, 33.6);
$hsp = $hit->next_hsp;
ok($hit->hsp->bits,$hsp->bits);

ok($hsp->get_aln->isa('Bio::Align::AlignI'));
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

my $out = new Bio::SearchIO(-writer => $writer,
			 -file   => ">searchio.out");
$out->write_result($result, 1);
ok(-e 'searchio.out');
my $writerhtml = new Bio::SearchIO::Writer::HTMLResultWriter();
my $outhtml = new Bio::SearchIO(-writer => $writerhtml,
				-file   => ">searchio.html");
$outhtml->write_result($result, 1);
ok(-e "searchio.html");

unlink 'searchio.out';
unlink 'searchio.html';

#test all the database accession number formats
$searchio = new Bio::SearchIO(-format => 'blast',
				 -file   => File::Spec->catfile(qw(t data testdbaccnums.out)) );
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
ok($hit->name, shift @$d);
ok($hit->accession, shift @$d);
my @accnums = $hit->each_accession_number;
foreach my $a (@accnums) {
	ok($a, shift @$d);
}
$d = shift @valid;
$hit = $result->next_hit;
ok($hit->name, shift @$d);
ok($hit->accession, shift @$d);
ok($hit->locus, shift @$d);

$hits_left = 23;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->accession, shift @$d);
    $hits_left--;
}
ok $hits_left, 0;

# Parse MEGABLAST

# parse the BLAST-like output
my $infile = Bio::Root::IO->catfile(qw(t data 503384.MEGABLAST.2));
my $in = new Bio::SearchIO(-file => $infile,
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

ok($r->query_name, '503384');
ok($r->query_description, '11337 bp 2 contigs');
ok($r->query_length, 11337);
ok($r->database_name, 'cneoA.nt ');
ok($r->database_letters, 17206226);
ok($r->database_entries, 4935);

$hits_left = 4;
while( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->raw_score, shift @$d);
    ok($hit->significance, shift @$d);
    
    my $hsp = $hit->next_hsp;
    ok($hsp->bits, shift @$d);
    ok($hsp->evalue, shift @$d);
    ok($hsp->score, shift @$d);
    ok($hsp->num_identical, shift @$d);
    ok($hsp->gaps('total'), shift @$d);
    ok($hsp->query->start, shift @$d);
    ok($hsp->query->end, shift @$d);
    ok($hsp->query->strand, shift @$d);
    ok($hsp->hit->start, shift @$d);
    ok($hsp->hit->end, shift @$d);
    ok($hsp->hit->strand, shift @$d);
    $hits_left--;
}
ok $hits_left, 0;

# parse another megablast format

$infile =  Bio::Root::IO->catfile(qw(t data 503384.MEGABLAST.0));

# this is megablast output type 0 
$in = new Bio::SearchIO(-file          => $infile,
			-report_format => 0,
			-format        => 'megablast'); 
$r = $in->next_result;
@dcompare = ( 
	      ['Contig634', 7620, 7941, 1, 1, 321, -1],
	      ['Contig1853', 6406, 6620, 1, 1691, 1905, 1],  
	      ['Contig3700', 8723,9434, 1, 4083, 4794, -1],
	      ['Contig3997', 1282, 1704, 1, 1546, 1968,-1 ],
	      );

ok($r->query_name, '503384');

while( my $hit = $r->next_hit ) {
    my $d = shift @dcompare;
    ok($hit->name, shift @$d);
    my $hsp = $hit->next_hsp;
    ok($hsp->query->start, shift @$d);
    ok($hsp->query->end, shift @$d);
    ok($hsp->query->strand, shift @$d);
    ok($hsp->hit->start, shift @$d);
    ok($hsp->hit->end, shift @$d);
    ok($hsp->hit->strand, shift @$d);
}
ok @dcompare, 0;


# Let's test RPS-BLAST

my $parser = new Bio::SearchIO(-format => 'blast',
			       -file   => Bio::Root::IO->catfile(qw(t data ecoli_domains.rpsblast)));

$r = $parser->next_result;
ok($r->query_name, 'gi|1786183|gb|AAC73113.1|');
ok($r->num_hits, 7);
$hit = $r->next_hit;
ok($hit->name, 'gnl|CDD|3919');
ok($hit->significance, 0.064);
ok($hit->raw_score, 28);
$hsp = $hit->next_hsp;
ok($hsp->query->start, 599);
ok($hsp->query->end,655);
ok($hsp->hit->start,23);
ok($hsp->hit->end,76);


# Test PSI-BLAST parsing

$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','psiblastreport.out'));

$result = $searchio->next_result;

ok($result->database_name, '/home/peter/blast/data/swissprot.pr');
ok($result->database_entries, 88780);
ok($result->database_letters, 31984247);

ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, qr/^2\.0\.14/);
ok($result->query_name, 'CYS1_DICDI');
ok($result->query_length, 343);
ok($result->get_statistic('kappa') == 0.0491);
ok($result->get_statistic('lambda') == 0.270);
ok($result->get_statistic('entropy') == 0.230);
ok($result->get_statistic('dbletters'), 31984247);
ok($result->get_statistic('dbentries'), 88780);
ok($result->get_statistic('effective_hsplength'), 49);
ok($result->get_statistic('effectivespace'), 8124403938);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_parameter('gapopen'), 11);
ok($result->get_parameter('gapext'), 1);

my @valid_hit_data = ( [ 'sp|P04988|CYS1_DICDI', 343, 'P04988', '0', 721],
		       [ 'sp|P43295|A494_ARATH', 313, 'P43295', '1e-75', 281],
		       [ 'sp|P25804|CYSP_PEA', 363, 'P25804', '1e-74', 278]);
my @valid_iter_data = ( [ 127, 127, 0, 109, 18, 0, 0, 0, 0],
			[ 157, 40, 117, 2, 38, 0, 109, 3, 5]);
my $iter_count = 0;
while( $iter = $result->next_iteration ) {
    $iter_count++;
    my $di = shift @valid_iter_data;
    ok($iter->number, $iter_count);

    ok($iter->num_hits, shift @$di);
    ok($iter->num_hits_new, shift @$di);
    ok($iter->num_hits_old, shift @$di);
    ok(scalar($iter->newhits_below_threshold), shift @$di);
    ok(scalar($iter->newhits_not_below_threshold), shift @$di);
    ok(scalar($iter->newhits_unclassified), shift @$di);
    ok(scalar($iter->oldhits_below_threshold), shift @$di);
    ok(scalar($iter->oldhits_newly_below_threshold), shift @$di);
    ok(scalar($iter->oldhits_not_below_threshold), shift @$di);

    my $hit_count = 0;
    if ($iter_count == 1) {
        while( $hit = $result->next_hit ) {
            my $d = shift @valid_hit_data;
            
            ok($hit->name, shift @$d);
            ok($hit->length, shift @$d);
            ok($hit->accession, shift @$d);
            ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
            ok($hit->bits, shift @$d );
            
            if( $hit_count == 1 ) {
                my $hsps_left = 1;
                while( my $hsp = $hit->next_hsp ){
                    ok($hsp->query->start, 32);
                    ok($hsp->query->end, 340);
                    ok($hsp->hit->start, 3);
                    ok($hsp->hit->end, 307);
                    ok($hsp->length('hsp'), 316);
                    ok($hsp->start('hit'), $hsp->hit->start);
                    ok($hsp->end('query'), $hsp->query->end);
                    ok($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
                    ok($hsp->evalue == 1e-75);
                    ok($hsp->score, 712);
                    ok($hsp->bits, 281);
                    ok(sprintf("%.1f",$hsp->percent_identity), 46.5);
                    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.4757);
                    ok(sprintf("%.3f",$hsp->frac_identical('hit')), 0.482);
                    ok($hsp->gaps, 18);
                    $hsps_left--;
                }
                ok $hsps_left, 0;
            }
            last if( $hit_count++ > @valid_hit_data );
        }
    }
}
ok @valid_hit_data, 0;
ok @valid_iter_data, 0;

# Test filtering

$searchio = new Bio::SearchIO ( '-format' => 'blast', 
                                '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'),
                                '-signif' => 1e-100);

@valid = qw(gb|AAC73113.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    ok($hit->name, shift @valid);
}

$searchio = new Bio::SearchIO ( '-format' => 'blast', 
                                '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'),
                                '-score' => 100);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1| gb|AAC76994.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    ok($hit->name, shift @valid);
}
ok @valid, 0;

$searchio = new Bio::SearchIO ( '-format' => 'blast', 
                                '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'),
                                '-bits' => 200);

@valid = qw(gb|AAC73113.1| gb|AAC76922.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    ok($hit->name, shift @valid);
}
ok @valid, 0;


my $filt_func = sub{ my $hit=shift; 
                     $hit->frac_identical('query') >= 0.31
                     };

$searchio = new Bio::SearchIO ( '-format' => 'blast', 
                                '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'),
                                '-hit_filter' => $filt_func);

@valid = qw(gb|AAC73113.1| gb|AAC76994.1|);
$r = $searchio->next_result;

while( my $hit = $r->next_hit ) {
    ok($hit->name, shift @valid);
}
ok @valid, 0;




# bl2seq parsing testing

# this is blastp bl2seq
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile(qw(t data
								   bl2seq.out)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, '');
ok($result->algorithm, 'BLASTP');
$hit = $result->next_hit;
ok($hit->name, 'ALEU_HORVU');
ok($hit->length, 362);
$hsp = $hit->next_hsp;
ok($hsp->score, 481);
ok($hsp->bits, 191);
ok(int $hsp->percent_identity, 34);
ok($hsp->evalue, '2e-53');
ok(int($hsp->frac_conserved*$hsp->length), 167);
ok($hsp->query->start, 28);
ok($hsp->query->end, 343);
ok($hsp->hit->start, 60);
ok($hsp->hit->end,360);
ok($hsp->gaps, 27);

# this is blastn bl2seq 
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data bl2seq.blastn.rev)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, '');
ok($result->algorithm, 'BLASTN');
ok($result->query_length, 180);
$hit = $result->next_hit;
ok($hit->length, 179);
ok($hit->name, 'human');
$hsp = $hit->next_hsp;
ok($hsp->score, 27);
ok($hsp->bits, '54.0');
ok(int $hsp->percent_identity, 88);
ok($hsp->evalue, '2e-12');
ok(int($hsp->frac_conserved*$hsp->length), 83);
ok($hsp->query->start, 94);
ok($hsp->query->end, 180);
ok($hsp->query->strand, 1);
ok($hsp->hit->strand, -1);
ok($hsp->hit->start, 1);
ok($hsp->hit->end,94);
ok($hsp->gaps, 7);

# this is blastn bl2seq 
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data bl2seq.blastn)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, '');
ok($result->query_length, 180);
ok($result->algorithm, 'BLASTN');
$hit = $result->next_hit;
ok($hit->name, 'human');
ok($hit->length, 179);
$hsp = $hit->next_hsp;
ok($hsp->score, 27);
ok($hsp->bits, '54.0');
ok(int $hsp->percent_identity, 88);
ok($hsp->evalue, '2e-12');
ok(int($hsp->frac_conserved*$hsp->length), 83);
ok($hsp->query->start, 94);
ok($hsp->query->end, 180);
ok($hsp->query->strand, 1);
ok($hsp->hit->strand, 1);
ok($hsp->hit->start, 86);
ok($hsp->hit->end,179);
ok($hsp->gaps, 7);

# this is blastp bl2seq
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data bl2seq.bug940.out)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, 'zinc');
ok($result->algorithm, 'BLASTP');
ok($result->query_description, 'finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0');
ok($result->query_length, 469);
$hit = $result->next_hit;
ok($hit->name, 'gi|4507985|');
ok($hit->description,'zinc finger protein 135 (clone pHZ-17) [Homo sapiens]. neo_id RS.ctg14243-000000.6.0');
ok($hit->length, 469);
$hsp = $hit->next_hsp;
ok($hsp->score, 1626);
ok($hsp->bits, 637);
ok(int $hsp->percent_identity, 66);
ok($hsp->evalue, '0.0');
ok(int($hsp->frac_conserved*$hsp->length), 330);
ok($hsp->query->start, 121);
ok($hsp->query->end, 469);
ok($hsp->hit->start, 1);
ok($hsp->hit->end,469);
ok($hsp->gaps, 120);
ok($hit->next_hsp); # there is more than one HSP here, 
                    # make sure it is parsed at least

# cannot distinguish between blastx and tblastn reports
# so we're only testing a blastx report for now

# this is blastx bl2seq
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data bl2seq.blastx.out)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, 'AE000111.1');
ok($result->query_description, 'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome');
ok($result->algorithm, 'BLASTX');
ok($result->query_length, 720);
$hit = $result->next_hit;
ok($hit->name, 'AK1H_ECOLI');
ok($hit->description,'P00561 Bifunctional aspartokinase/homoserine dehydrogenase I (AKI-HDI) [Includes: Aspartokinase I ; Homoserine dehydrogenase I ]');
ok($hit->length, 820);
$hsp = $hit->next_hsp;
ok($hsp->score, 634);
ok($hsp->bits, 248);
ok(int $hsp->percent_identity, 100);
ok($hsp->evalue, '2e-70');
ok(int($hsp->frac_conserved*$hsp->length), 128);
ok($hsp->query->start, 1);
ok($hsp->query->end, 384);
ok($hsp->hit->start, 1);
ok($hsp->hit->end,128);
ok($hsp->gaps, 0);
ok($hsp->query->frame,0);
ok($hsp->hit->frame,0);
ok($hsp->query->strand,-1);
ok($hsp->hit->strand,0);

# this is tblastx bl2seq (self against self)
$searchio = new Bio::SearchIO(-format => 'blast',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data bl2seq.tblastx.out)));
$result = $searchio->next_result;
ok($result);
ok($result->query_name, 'Escherichia');
ok($result->algorithm, 'TBLASTX');
ok($result->query_description, 'coli K-12 MG1655 section 1 of 400 of the complete genome');
ok($result->query_length, 720);
$hit = $result->next_hit;
ok($hit->name, 'gi|1786181|gb|AE000111.1|AE000111');

ok($hit->description,'Escherichia coli K-12 MG1655 section 1 of 400 of the complete genome');
ok($hit->length, 720);
$hsp = $hit->next_hsp;
ok($hsp->score, 1118);
ok($hsp->bits, 515);
ok(int $hsp->percent_identity, 95);
ok($hsp->evalue, '1e-151');
ok(int($hsp->frac_conserved*$hsp->length), 229);
ok($hsp->query->start, 1);
ok($hsp->query->end, 720);
ok($hsp->hit->start, 1);
ok($hsp->hit->end,720);
ok($hsp->gaps, 0);
ok($hsp->query->frame,0);
ok($hsp->hit->frame,0);
ok($hsp->query->strand,1);
ok($hsp->hit->strand,1);

# this is NCBI tblastn
$searchio = new Bio::SearchIO(-format => 'blast',
										-file   => Bio::Root::IO->catfile
										(qw(t data tblastn.out)));
$result = $searchio->next_result;
ok($result);
ok($result->algorithm, 'TBLASTN');
$hit = $result->next_hit;
ok($hit->name,'gi|10040111|emb|AL390796.6|AL390796');

# test blasttable output
my @eqset = qw( 
		
		c200-vs-yeast.BLASTN.m9);
$searchio = new Bio::SearchIO(-file => Bio::Root::IO->catfile
			      (qw(t data c200-vs-yeast.BLASTN)),
			      -format => 'blast');
$result = $searchio->next_result;
ok($result);
my %ref = &result2hash($result);
ok( scalar keys %ref, 67);
$searchio = new Bio::SearchIO(-file => Bio::Root::IO->catfile
			      (qw(t data c200-vs-yeast.BLASTN.m8)),
			      -program_name => 'BLASTN',
			      -format => 'blasttable');
$result = $searchio->next_result;
my %tester = &result2hash($result);
ok( scalar keys %tester, 67);
foreach my $key ( sort keys %ref ) {
    ok($tester{$key}, $ref{$key},$key);
}      




# Test Blast parsing with B=0 (WU-BLAST)
$searchio = new Bio::SearchIO(-file   => Bio::Root::IO->catfile
			      (qw(t data no_hsps.blastp)),
			      -format => 'blast');
$result = $searchio->next_result;
ok($result->query_name, 'mgri:MG00189.3');
$hit = $result->next_hit;
ok($hit->name, 'mgri:MG00189.3');
ok($hit->description, 'hypothetical protein 6892 8867 +');
ok($hit->score, 3098);
ok($hit->significance, '0.');

$hit = $result->next_hit;
ok($hit->name, 'fgram:FG01141.1');
ok($hit->description, 'hypothetical protein 47007 48803 -');
ok($hit->score, 2182);
ok($hit->significance, '4.2e-226');
ok($result->num_hits, 415);
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
    ok(Bio::SearchIO->_guess_format($file),$expformat, "$expformat for $file");
}


# Test Wes Barris's reported bug when parsing blastcl3 output which
# has integer overflow

$searchio = new Bio::SearchIO(-file => Bio::Root::IO->catfile
			      (qw(t data hsinsulin.blastcl3.blastn)),
			      -format => 'blast');
$result = $searchio->next_result;
ok($result->query_name, 'human');
ok($result->database_letters(), '-24016349'); 
# this is of course not the right length, but is the what blastcl3 
# reports, the correct value is
ok($result->get_statistic('dbletters'),'192913178');
ok($result->get_statistic('dbentries'),'1867771');


# test for links and groups being parsed out of WU-BLAST properly
$searchio = Bio::SearchIO->new(-format => 'blast',
			       -file   => Bio::Root::IO->catfile
			       (qw(t data brassica_ATH.WUBLASTN) ));
ok($result = $searchio->next_result);
ok($hit = $result->next_hit);
ok($hsp = $hit->next_hsp);
ok($hsp->links,'(1)-3-2');
ok($hsp->query->strand, 1);
ok($hsp->hit->strand, 1);
ok($hsp->hsp_group, '1');

## Web blast result parsing

$searchio = Bio::SearchIO->new(-format => 'blast',
			       -file   => Bio::Root::IO->catfile
			       (qw(t data catalase-webblast.BLASTP)));
ok($result = $searchio->next_result);
ok($hit = $result->next_hit);
ok($hit->name, 'gi|40747822|gb|EAA66978.1|', 'full hit name');
ok($hit->accession, 'EAA66978', 'hit accession');
ok($hsp = $hit->next_hsp);
ok($hsp->query->start, 1, 'query start');
ok($hsp->query->end, 528, 'query start');

# tests for new BLAST 2.2.13 output
$searchio = new Bio::SearchIO(-format => 'blast',
							  -file   => Bio::Root::IO->catfile
							  (qw(t data new_blastn.txt)));

$result = $searchio->next_result;
ok($result->database_name, 'All GenBank+EMBL+DDBJ+PDB sequences (but no EST, STS,GSS,environmental samples or phase 0, 1 or 2 HTGS sequences)');
ok($result->database_entries, 3742891);
ok($result->database_letters, 16670205594);
ok($result->algorithm, 'BLASTN');
ok($result->algorithm_version, '2.2.13 [Nov-27-2005]');
ok($result->query_name, 'pyrR,');
ok($result->query_length, 558);
ok($result->get_statistic('kappa'), '0.711');
ok($result->get_statistic('kappa_gapped'), '0.711');
ok($result->get_statistic('lambda'), '1.37');
ok($result->get_statistic('lambda_gapped'), '1.37');
ok($result->get_statistic('entropy'), '1.31');
ok($result->get_statistic('entropy_gapped'), '1.31');
ok($result->get_statistic('dbletters'), '-509663586');
ok($result->get_statistic('dbentries'), 3742891);
ok($result->get_statistic('effective_hsplength'), undef);
ok($result->get_statistic('effectivespace'), undef);
ok($result->get_parameter('matrix'), 'blastn matrix:1 -3');
ok($result->get_parameter('gapopen'), 5);
ok($result->get_parameter('gapext'), 2);
ok($result->get_statistic('S2'), '60');
ok($result->get_statistic('S2_bits'), '119.4');
ok($result->get_parameter('expect'), '1e-23');
ok($result->get_statistic('num_extensions'), '117843');


@valid = ( [ 'gi|41400296|gb|AE016958.1|', 4829781, 'AE016958', '6e-059', 236],
	      [ 'gi|54013472|dbj|AP006618.1|', 6021225, 'AP006618', '4e-026', 127],
	      [ 'gi|57546753|dbj|BA000030.2|', 9025608, 'BA000030', '1e-023', 119]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok(sprintf("%g",$hit->significance), sprintf("%g",shift @$d) );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
        my $hsps_left = 1;
        while( my $hsp = $hit->next_hsp ) {
            ok($hsp->query->start, 262);
            ok($hsp->query->end, 552);
            ok($hsp->hit->start, 1166897);
            ok($hsp->hit->end, 1167187);
            ok($hsp->length('hsp'), 291);
            ok($hsp->start('hit'), $hsp->hit->start);
            ok($hsp->end('query'), $hsp->query->end);
            ok($hsp->strand('sbjct'), $hsp->subject->strand);# alias for hit
            ok($hsp->evalue == 6e-59);
            ok($hsp->score, 119);
            ok($hsp->bits,236);	    	    
            ok(sprintf("%.2f",$hsp->percent_identity), 85.22);
            ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.8522);
            ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.8522);
            ok($hsp->gaps, 0);
            $hsps_left--;
        }
        ok $hsps_left, 0;
    }
    last if( $count++ > @valid );
}
ok @valid, 0;

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


