# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;

use vars qw($SKIPXML $LASTXMLTEST); 
use strict;
use lib '.';

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use vars qw($NTESTS);
    $NTESTS = 248;
    $LASTXMLTEST = 49;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 

    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	$SKIPXML = 1;
	print STDERR "XML::Parser::PerlSAX not loaded. This means SearchIO::blastxml test cannot be executed. Skipping\n";
	foreach ( 1..$LASTXMLTEST ) {
	    skip('No XML::Parser::PerlSAX loaded',1);
	}
    }
}

if( $error == 1 ) {
    exit(0);
}

use Bio::SearchIO;
use Bio::Root::IO;

ok(1);

my ($searchio, $result,$hit,$hsp);
if( ! $SKIPXML ) {
    # test with RPSBLAST data first 
    $searchio = new Bio::SearchIO ('-format' => 'blastxml',
				   '-file'   => Bio::Root::IO->catfile('t','data','ecoli_domains.rps.xml'));
    
    $result = $searchio->next_result;
    ok($result);    
    ok($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam');
    ok($result->query_name,'gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide [Escherichia coli]');
    ok($result->query_length, 21);
    ok($result->algorithm, 'BLASTP');
    ok($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');

    ok($result->available_parameters, 8);
    ok($result->get_parameter('gapext'), 1);
    ok($result->available_statistics, 7);
    ok($result->get_statistic('lambda'), 0.267);

# this result actually has a hit
    $result = $searchio->next_result;
    $hit = $result->next_hit;
    ok($hit->name, 'gnl|Pfam|pfam00742');
    ok($hit->description(), 'HomoS_dh, HomoS dehydrogenase');
    ok($hit->accession, 'pfam00742');
    ok($hit->length, 310);

    $hsp = $hit->next_hsp;
    ok($hsp->p, 1.46134e-90);
    ok($hsp->evalue, 1.46134e-90);
    ok($hsp->score, 838);
    ok($hsp->bits,327.405);
    ok($hsp->query->start, 498);
    ok($hsp->query->end,815);
    ok($hsp->hit->start, 3);
    ok($hsp->hit->end, 310);
    ok($hsp->query->frame,0);
    ok($hsp->hit->frame,0);
    ok(sprintf("%.2f", $hsp->percent_identity), 37.73);
    ok(sprintf("%.4f", $hsp->frac_identical('hit')), 0.3968);
    ok($hsp->frac_identical('query'), 0.15);
    ok($hsp->query->frac_identical, 0.15);

    ok(sprintf("%.4f",$hsp->frac_conserved('total')),0.3773);
    ok(sprintf("%.4f",$hsp->frac_conserved('hit')),0.3968);
    ok($hsp->frac_conserved('query'),0.15);
    ok($hsp->gaps('total'), 26);
    ok($hsp->length('hsp'), 326);
    ok($hsp->query_string, 'LRVCGVANSKALLTNVHGLNLENWQEELAQAKEPF-NLGRLIRLVKEYHLLN----PVIVDCTSSQAVAD-QYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE-GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARET-GRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS');
    ok($hsp->hit_string, 'GVVTGITDSREMLLSRIGLPLEIWKVALRDLEKPRKDLGKLDLTDDAFAVVDDPDIDVVVELTGGIEVARELYLDALEEGKHVVTANKALNASHGDEYLAL---AEKSGVDVLYEAAVAGGIPIIKTLRELLATGDRILKIEGIFNGTTNFILSEMDEKGLPFSDVLAEAQELGYTEADPRDDVEGIDAARKLAILARIAFGIELELDDVYVEGISPITAEDISSADEFGYTLKLLDEAMRQRVEDAESGGEVLRYPTLIPE-------------DHPLASVKGSDNAVAVEGEAYG--PLMFYGPGAGAEPTASAVVADIVRIAR');
    ok($hsp->homology_string, '  V G+ +S+ +L +  GL LE W+  L   ++P  +LG+L      + +++     V+V+ T    VA   Y D L EG HVVT NK  N S  D Y  L   AEKS    LY+  V  G+P+I+ L+ LL  GD ++K  GI +G+ ++I  ++DE G+ FS+    A+E+GYTE DPRDD+ G+D ARKL ILAR   G ELEL D+ +E + P           F   L  LD+    RV  A   G+VLRY   I E             + PL  VK  +NA+A     Y   PL+  G GAG + TA+ V AD++R   ');

# one more 
    $hit = $result->next_hit;
    ok($hit);

    while( $result = $searchio->next_result ) { ok($result); }


    $searchio = new Bio::SearchIO(-format => 'blastxml', -file => Bio::Root::IO->catfile('t','data','plague_yeast.bls.xml'));

    $result = $searchio->next_result;

    ok($result->database_name, 'yeast.aa');
    ok($result->query_name, 'gi|5763811|emb|CAB53164.1| putative transposase [Yersinia pestis]');
    ok($result->query_length, 340);

    $hit = $result->next_hit;
    ok(! $hit);

}

$searchio = new Bio::SearchIO ('-format' => 'blast',
				  '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'));

$result = $searchio->next_result;

ok($result->database_name, 'ecoli.aa');
ok($result->database_entries, 1358990);
ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, '2.1.3');
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_length, 820);
ok($result->get_statistic('kappa')== 0.041);
ok($result->get_statistic('lambda'), 0.267);
ok($result->get_statistic('entropy') == 0.14);
ok($result->get_statistic('dblength'), 1358990);
ok($result->get_statistic('dbnum'), 4289);
ok($result->get_statistic('hsplength'), 47);
ok($result->get_statistic('effectivespace'), 894675611);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_parameter('gapopen'), 11);
ok($result->get_parameter('gapext'), 1);

my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113'],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922'],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994']);
my $count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, $d->[0]);
    ok($hit->length, $d->[1]);
    ok($hit->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 820);
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 820);
	    ok($hsp->length('hsp'), 820);
	    ok($hsp->p == 0.0);
	    ok($hsp->evalue == 0.0);
	    ok($hsp->score, 4058);
	    ok($hsp->bits,1567);	    	    
	    ok(sprintf("%.2f",$hsp->percent_identity), 98.29);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.9829);
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.9829);
	    ok($hsp->gaps, 0);	    
	}
    }
    last if( $count++ > @valid );
}

$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.wublastn'));

$result = $searchio->next_result;

ok($result->database_name, 'ecoli.aa');
ok($result->database_entries, 1358990);
ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, '2.0MP-WashU');
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_length, 820);
ok($result->get_statistic('kappa'), 0.136);
ok($result->get_statistic('lambda'), 0.319);
ok($result->get_statistic('entropy'), 0.384);
ok($result->get_statistic('dblength'), 1358990);
ok($result->get_statistic('dbnum'), 4289);
ok($result->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113'],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922'],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994']);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, $d->[0]);
    ok($hit->length, $d->[1]);
    ok($hit->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 820);
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 820);
	    ok($hsp->length('hsp'), 820);
	    ok($hsp->p == 0.0);
	    ok($hsp->evalue == 0.0);
	    ok($hsp->score, 4141);
	    ok($hsp->bits,1462.8);	    	    
	    ok($hsp->percent_identity, 100);
	    ok($hsp->frac_identical('query'), 1.00);
	    ok($hsp->frac_identical('hit'), 1.00);
	    ok($hsp->gaps, 0);	    
	}
    }
    last if( $count++ > @valid );
}

# test tblastx 
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

$result = $searchio->next_result;
ok($result->database_name, 'ecoli.nt');
ok($result->database_entries, 4662239);
ok($result->algorithm, 'TBLASTX');
ok($result->algorithm_version, '2.1.2');
ok($result->query_name, qr/HUMBETGLOA Human haplotype C4 beta-globin gene, complete cds./);
ok($result->query_length, 3002);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('lambda'), 0.318);
ok($result->get_statistic('entropy'), 0.401);
ok($result->get_statistic('dblength'), 4662239);
ok($result->get_statistic('dbnum'), 400);
ok($result->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AE000479.1|AE000479', 10934, 'AE000479'],
	   [ 'gb|AE000302.1|AE000302', 10264, 'AE000302'],
	   [ 'gb|AE000277.1|AE000277', 11653, 'AE000277']);
$count = 0;

while( $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, $d->[0]);
    ok($hit->length, $d->[1]);
    ok($hit->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1057);
	    ok($hsp->query->end, 1134);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->end, 5893);
	    ok($hsp->hit->start, 5816);
	    ok($hsp->hit->strand, -1);
	    ok($hsp->length('hsp'), 26);
	    ok($hsp->p == 0.13);
	    ok($hsp->evalue == 0.13);
	    ok($hsp->score, 67);
	    ok($hsp->bits,33.6);
	    ok(sprintf("%.2f",$hsp->percent_identity), 42.31);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.0037);
	    ok(sprintf("%.3f",$hsp->frac_identical('hit')), 0.001);
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 1);
	    ok($hsp->gaps, 0);	    
	    ok($hsp->query_string, 'SAYWSIFPPLGCWWSTLGPRGSLSPL');
	    ok($hsp->hit_string, 'AAVWALFPPVGSQWGCLASQWRTSPL');
	    ok($hsp->homology_string, '+A W++FPP+G  W  L  +   SPL');
	}
    }
    last if( $count++ > @valid );
}

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/HUMBETGLOA.FASTA');
$result = $searchio->next_result;
ok($result->database_name, qr/dros_clones.2.5/);
ok($result->database_entries, 112936249);
ok($result->algorithm, 'FASTA');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, qr/HUMBETGLOA Human haplotype C4 beta-globin gene, complete cds./);
ok($result->query_length, 3002);
ok($result->get_parameter('gapopen'), -16);
ok($result->get_parameter('gapext'), -4);
ok($result->get_parameter('ktup'), 6);

ok($result->get_statistic('lambda'), 0.0823);
ok($result->get_statistic('dblength'), 112936249);
ok($result->get_statistic('dbnum'), 657);

@valid = ( [ 'BACR21I23', 73982, 'BACR21I23'],
	   [ 'BACR40P19', 73982, 'BACR40P19'],
	   [ 'BACR30L17', 32481, 'BACR30L17']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, $d->[0]);
    ok($hit->length, $d->[1]);
    ok($hit->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 31);
	    ok($hsp->query->end, 289);
	    ok($hsp->query->strand, -1);
	    ok($hsp->hit->end, 65167);
	    ok($hsp->hit->start, 64902);
	    ok($hsp->hit->strand, 1);
	    ok($hsp->length('hsp'), 267);
	    ok($hsp->p == 0.017);
	    ok($hsp->evalue == 0.017);
	    ok($hsp->score, 134.5);
	    ok($hsp->bits,44.2);
	    ok(sprintf("%.1f",$hsp->percent_identity), 57.3);
	    ok(sprintf("%.3f",$hsp->frac_identical('query')), 0.051); 
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.0021); 
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 0);
	    ok($hsp->gaps, 159);
	    ok($hsp->gaps('query'), 1);
	    ok($hsp->gaps('hit'),1);
	    ok($hsp->query_string, 'GATTAAAACCTTCTGGTAAGAAAAGAAAAAATATATATATATATATATGTGTATATGTACACACATACATATACATATATATGCATTCATTTGTTGTTGTTTTTCTTAATTTGCTCATGCATGCTA----ATAAATTATGTCTAAAAATAGAAT---AAATACAAATCAATGTGCTCTGTGCATTA-GTTACTTATTAGGTTTTGGGAAACAAGAGGTAAAAAACTAGAGACCTCTTAATGCAGTCAAAAATACAAATAAATAAAAAGTCACTTACAACCCAAAGTGTGACTATCAATGGGGTAATCAGTGGTGTCAAATAGGAGGT');
	    ok($hsp->hit_string, 'GATGTCCTTGGTGGATTATGGTGTTAGGGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATATAAAATATAATATAAAATATAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATAT-AATATAAAATATAAAATAAAATATAATATAAAATATAATATAAAATATAATATAAAATATAATATAAAATA');
	    ok($hsp->homology_string, '                              :::::::::::::::::: : ::::: :: : : ::: ::::: ::::::::  ::  :: : :   : : : : :  ::    : :: ::   ::    : ::: :::     :::::: :::   ::::: ::  :::  :    :    : ::   :::  : ::   : :   : : :: :   :: : : :: : :       ::  : : ::: ::: ::  ::::: ::: : :  :: ::   ::: : : : ::: ::                                                               ');
	}
    }
    last if( $count++ > @valid );
} 

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/cysprot1.FASTA');
$result = $searchio->next_result;
ok($result->database_name, qr/ecoli.aa/);
ok($result->database_entries, 1358987);
ok($result->algorithm, 'FASTA');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, 'CYS1_DICDI');
ok($result->query_length, 343);
ok($result->get_parameter('gapopen'), -12);
ok($result->get_parameter('gapext'), -2);
ok($result->get_parameter('ktup'), 2);

ok($result->get_statistic('lambda'), 0.1456);
ok($result->get_statistic('dblength'), 1358987);
ok($result->get_statistic('dbnum'), 4289);


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309'],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148'],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494']);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, $d->[0]);
    ok($hit->length, $d->[1]);
    ok($hit->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 125);
	    ok($hsp->query->end, 305);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->start, 2);
	    ok($hsp->hit->end, 181);
	    ok($hsp->hit->strand, 1);
	    ok($hsp->length('hsp'), 188);
	    ok($hsp->p == 1.2);
	    ok($hsp->evalue == 1.2);
	    ok($hsp->score, 109.2);
	    ok($hsp->bits,29.2);
	    ok(sprintf("%.2f",$hsp->percent_identity), 23.94);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.1312);
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.0879);
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 0);
	    ok($hsp->gaps, 49);	    
	    ok($hsp->query_string, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
	    ok($hsp->hit_string, 'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
	    ok($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                                                                                  ');
	}
    }
    last if( $count++ > @valid );
} 


# TODO: Flesh this test out!
$searchio = new Bio::SearchIO ('-format' => 'psiblast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

 $result = $searchio->next_result;

 ok($result);

