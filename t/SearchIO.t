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
    $NTESTS = 601;
    $LASTXMLTEST = 50;
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
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;

ok(1);
my ($searchio, $result,$hit,$hsp);
if( ! $SKIPXML ) {
    # test with RPSBLAST data first 
    $searchio = new Bio::SearchIO ('-tempfile' => 1,
				   '-format' => 'blastxml',
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
    ok($hsp->pvalue, undef);
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
    
    while( $result = $searchio->next_result ) { ok($result); }


    $searchio = new Bio::SearchIO(-format => 'blastxml', 
				  -file => Bio::Root::IO->catfile('t','data','plague_yeast.bls.xml'));

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
ok($result->database_entries, 4289);
ok($result->database_letters, 1358990);

ok($result->algorithm, 'BLASTP');
ok($result->algorithm_version, qr/^2\.1\.3/);
ok($result->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($result->query_length, 820);
ok($result->get_statistic('kappa')== 0.041);
ok($result->get_statistic('lambda'), 0.267);
ok($result->get_statistic('entropy') == 0.14);
ok($result->get_statistic('dbletters'), 1358990);
ok($result->get_statistic('dbentries'), 4289);
ok($result->get_statistic('hsplength'), 47);
ok($result->get_statistic('effectivespace'), 894675611);
ok($result->get_parameter('matrix'), 'BLOSUM62');
ok($result->get_parameter('gapopen'), 11);
ok($result->get_parameter('gapext'), 1);

my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113.1', '0.0', 1567],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922.1', '1e-91', 332],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994.1', '3e-47', 184]);
my $count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 820);
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 820);
	    ok($hsp->length('hsp'), 820);
	    
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

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113.1', '0.0', 4141],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922.1', '3.1e-86', 844],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994.1', '2.8e-47', 483]);
$count = 0;
while( $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
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
	}
    }
    last if( $count++ > @valid );
}

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
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1057);
	    ok($hsp->query->end, 1134);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->end, 5893);
	    ok($hsp->hit->start, 5816);
	    ok($hsp->hit->strand, -1);
	    ok($hsp->length('hsp'), 26);
	    
	    ok($hsp->evalue == 0.13);
	    ok($hsp->score, 67);
	    ok($hsp->bits,33.6);
	    ok(sprintf("%.2f",$hsp->percent_identity), 42.31);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), '0.1410');
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.1410');
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 1);
	    ok($hsp->gaps, 0);	    
	    ok($hsp->query_string, 'SAYWSIFPPLGCWWSTLGPRGSLSPL');
	    ok($hsp->hit_string, 'AAVWALFPPVGSQWGCLASQWRTSPL');
	    ok($hsp->homology_string, '+A W++FPP+G  W  L  +   SPL');
	    ok(join(' ', $hsp->seq_inds('query', 'nomatch',1)), '355 364 365 367 368 370 371 373-375');
	}
    }
    last if( $count++ > @valid );
}

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/HUMBETGLOA.FASTA');
$result = $searchio->next_result;
ok($result->database_name, qr/dros_clones.2.5/);
ok($result->database_letters, 112936249);
ok($result->database_entries, 657);
ok($result->algorithm, 'FASTN');
ok($result->algorithm_version, '3.3t08');
ok($result->query_name, qr/HUMBETGLOA Human haplotype C4 beta-globin gene, complete cds./);
ok($result->query_length, 3002);
ok($result->get_parameter('gapopen'), -16);
ok($result->get_parameter('gapext'), -4);
ok($result->get_parameter('ktup'), 6);

ok($result->get_statistic('lambda'), 0.0823);
ok($result->get_statistic('dbletters'), 112936249);
ok($result->get_statistic('dbentries'), 657);

@valid = ( [ 'BACR21I23', 73982, 'BACR21I23', '0.017', 44],
	   [ 'BACR40P19', 73982, 'BACR40P19', '0.017', 44],
	   [ 'BACR30L17', 32481, 'BACR30L17', '0.018', 44]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );
    
    if( $count == 0 ) {
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
	}
    }
    last if( $count++ > @valid );
} 

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/cysprot1.FASTA');
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


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309', 1.2, 29],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148', 2.1, 27],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494',  2.1, 26]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
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
	    ok($hsp->hit_string, '                             MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
	    ok($hsp->homology_string, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                      '.' 'x60);
	}
    }
    last if( $count++ > @valid );
} 

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/cysprot_vs_gadfly.FASTA');
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
	     'FBan0006692', 3.1e-59, 228],
	   [ 'CG11459|FBgn0037396|pp-CT28891|FBan0011459', 336, 
	     'FBan0011459', 6.4e-41,  167],
	   [ 'CG4847|FBgn0034229|pp-CT15577|FBan0004847', 390, 
	     'FBan0004847',  2.5e-40, 165]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1373);
	    ok($hsp->query->end, 1706);
	    ok($hsp->query->strand, 0);
	    ok($hsp->hit->start, 5);
	    ok($hsp->hit->end, 341);
	    ok($hsp->hit->strand, 0);
	    ok($hsp->length('hsp'), 345);	    
	    ok($hsp->evalue == 3.1e-59 );
	    ok($hsp->score, 1170.6);
	    ok($hsp->bits,227.8);
	    ok(sprintf("%.2f",$hsp->percent_identity), 53.04);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.5479);
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.5430');
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 0);
	    ok($hsp->gaps('query'), 11);
	    ok($hsp->gaps, 194);
	    ok($hsp->query_string, 'SNWGNNGYFLIERGKNMCGLAACASYPIPQVMNPTLILAAFCLGIASATLTFDHSLEAQWTKWKAMHNRLY-GMNEEGWRRAVWEKNMKMIELHNQEYREGKHSFTMAMNAFGDMTSEEFRQVMNGFQ---NRKPR------KGKVFQEPLFYEAPRSVDWREKGYVTPVKNQGQCGSCWAFSATGALEGQMFRKTGRLISLSEQNLVDCSGPQGNEGCNGGLMDYAFQYVQDNGGLDSEESYPYEATEESCKYNPKYSVANDTGFVDIPK-QEKALMKAVATVGPISVAIDAGHESFLFYKEGIYFEPDCSSEDMDHGVLVVGYGFESTESDNNKYWLVKNSWGEEWGMGGYVKMAKDRRNHCGIASAASYPTVMTPLLLLAVLCLGTALATPKFDQTFNAQWHQWKSTHRRLYGTNEE');

	}
    }
    last if( $count++ > @valid );
} 


# test on TFASTXY
$searchio = new Bio::SearchIO(-format => 'fasta',
			      -file   => 't/data/5X_1895.FASTXY');
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


@valid = ( [ 'NR_SC:SW-YNN2_YEAST', 1056, 'NR_SC:SW-YNN2_YEAST','1.6e-154', 547],
	   [ 'NR_SC:SW-MPCP_YEAST', 311, 'NR_SC:SW-MPCP_YEAST', '1.3e-25', 117],
	   [ 'NR_SC:SW-YEO3_YEAST', 300, 'NR_SC:SW-YEO3_YEAST', '5.7e-05', 48]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;

    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 2180);
	    ok($hsp->query->end, 5623);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->start, 3);
	    ok($hsp->hit->end, 1053);
	    ok($hsp->hit->strand, 0);
	    ok($hsp->length('hsp'), 1165);

	    ok($hsp->evalue == 1.6e-154);
	    ok($hsp->score, 2877.6);
	    ok($hsp->bits,'547.0');
	    ok(sprintf("%.2f",$hsp->percent_identity), 51.67);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.1748);
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.5728);
	    ok($hsp->query->frame(), 0);
	    ok($hsp->hit->frame(), 0);
	    ok($hsp->gaps, 678);	    
	    ok($hsp->query_string, 'RKQLDPRIPALINNGVKANHRSFFVMVGDKGRDQVCPGMQAAMRFD*HRCR/LVNLHFLLSQARVSSRPSVLWCYKKD-LGFTT*VAASENLQQTIYFRPIATSHRKKREAKIKRDVKRGIRDANEQDPFELFVTVTDIRYTYYKDSAKILGQTFGMLVLQDYEAITPNLLARTIETVEGGGIVVLLLKTMSSLKQLYAMAM/DKL*CRDGVE*SDFS*LLI*DVHSRYRTDAHQFVQPRFNERFILSLGSNPDCLVLDDELNVLPLSKGKDIQIGKAGEEDDRGRKRKAEELKEMKENLEGVDIVGSLAKLAKTVDQAKAILTFVEAISEKNLSSTVALTAGRGRGKSAALGLAIGAALAHDYSNIFVTSPDPENLKTLFEFVFKALDALGYEEHIDYDVVQSTNPDFKKAIVRVNIFRGHRQTIQYISPEDSHVLGQAELVIIDEAAAIPLPLVRKLIGPYLVFMASTINGYEGTGRSLSIKLIQQLREQTRPSITKDSENAAASSAGSSSKAAAAGRSGAGLVRSLREIKLDEPIRYSPGDNVEKWLNNLLCLDATIVSK---SIQGCPHPSKCELYYVNRDTLFSYHPASEVFLQRMMALYVASHYKNSPNDLQMLSDAPAHHLFVLLPPIDEND-NTLPDPLVVLQVALEGNISREAILKEMAQSGMRSSGDMIPWIISTQFQDNDFATLSGARVVRIATHPDYARMGYGSRAMEALESFYNGTSYNFDDVPVDMGESFAD\VPRSDL*VTSFIPFPQNRTSTECVSQNANLQNDTIAIRDPSRMPPLLQRLSERKPETLDYLGVSFGLTRDLLRFWKKGGFTPLYASQKENALTGEYTFVMLKVLASAGGGGEWLGAFAQGMSCLLLQDEVHMGND*RL*TDFRQRFMNLLSYEAFKKFDASIALSILESTVPRNSPSPAP----KLLTNTELSSLLTPFDIKRLESYADSMLDYHVVLDLVPTIASLFFGKRLETS--LPPAQQAILLALGLQRKNVEALENELGITSTQTLALFGKVLRKMTKSLEDIRKASIASELP-----AEPTLAGRSANGSNKFVALQQTIEQDLADSAVQLNGEDDDASKKEQRELLNTLNMEEFAI-DQGGDWTEAEKQVERLASGKGGTRLSSTVSVKVDKLDD\AKRRRRRARMRVPRMRRR');
	    ok($hsp->hit_string, 'KKAIDSRIPSLIRNGVQTKQRSIFVIVGDRARNQ------------------LPNLHYLMMSADLKMNKSVLWAYKKKLLGFT--------------------SHRKKRENKIKKEIKRGTREVNEMDPFESFISNQNIRYVYYKESEKILGNTYGMCILQDFEALTPNLLARTIETVEGGGIVVILLKSMSSLKQLYTMTM-D--------------------VHARYRTEAHGDVVARFNERFILSLGSNPNCLVVDDELNVLPLSGAKNVKPLPPKEDDELPPKQL--ELQELKESLEDVQPAGSLVSLSKTVNQAHAILSFIDAISEKTLNFTVALTAGRGRGKSAALGISIAAAVSHGYSNIFVTSPSPENLKTLFEFIFKGFDALGYQEHIDYDIIQSTNPDFNKAIVRVDIKRDHRQTIQYIVPQDHQVLGQAELVVIDEAAAIPLPIVKNLLGPYLVFMASTINGYEGTGRSLSLKLIQQLRNQNNTSGRESTQTAVVSRDNKEKDSHLHSQS-----RQLREISLDEPIRYAPGDPIEKWLNKLLCLDVTLIKNPRFATRGTPHPSQCNLFVVNRDTLFSYHPVSENFLEKMMALYVSSHYKNSPNDLQLMSDAPAHKLFVLLPPIDPKDGGRIPDPLCVIQIALEGEISKESVRNSLSR-GQRAGGDLIPWLISQQFQDEEFASLSGARIVRIATNPEYASMGYGSRAIELLRDYFEGKF-------TDMSE---D-VRPKDYSI--------KRVSDKELAKT-NLLKDDVKLRDAKTLPPLLLKLSEQPPHYLHYLGVSYGLTQSLHKFWKNNSFVPVYLRQTANDLTGEHTCVMLNVLE--GRESNWLVEFAK---------------------DFRKRFLSLLSYD-FHKFTAVQALSVIESSKKAQDLSDDEKHDNKELTRTHLDDIFSPFDLKRLDSYSNNLLDYHVIGDMIPMLALLYFGDKMGDSVKLSSVQSAILLAIGLQRKNIDTIAKELNLPSNQTIAMFAKIMRKMSQYFRQLLSQSIEETLPNIKDDAIAEMDGEEIKNYNAAEALDQ-MEEDLEEAG----SEAVQAMREKQKELINSLNLDKYAINDNSEEWAESQKSLEIAAKAKGVVSLKTGKKRTTEKAED-IYRQEMKA-MKKPRKSKK');
	    ok($hsp->homology_string, '.: .: :::.:: :::....::.::.:::..:.:                  : :::.:. .: ..   :::: :::  ::::                    ::::::: :::...::: :..::.:::: :..  .:::.:::.: ::::.:.:: .:::.::.:::::::::::::::::::.:::.::::::::.:.: :                    ::.::::.::  :  ::::::::::::::.:::.:::::::::: .:...     :.:.   :.   ::.:.::.:: :. .:::..:.:::.::.:::.:..:::::.:. :::::::::::::::::..:.::..: :::::::::.::::::::::.::..:::::.::::::..:::::::.::::::.: : :::::::: :.: .::::::::.::::::::::.:..:.::::::::::::::::::::::.:::::::.:.  :  .....:..:  .. . .   ..:     :.::::.:::::::.::: .:::::.:::::.:....   . .: ::::.:.:. :::::::::::.:: ::..::::::.:::::::::::..::::::.::::::::: .: . .:::: :.:.::::.::.:.. . ... :.:..::.:::.:: ::::..::.:::::.:::::.:.:: :::::::.: :.....:         .::.:   : :  .:  .        .:.: . .... :: .: . .:: . .:::: .:::. :. : :::::.:::..: .:::...:.:.:  :  : ::::.: :::.::   :  ..::  ::.                     :::.::..::::. :.:: :  :::..::.   .. :       : :: :.:.....:::.:::.::....:::::. :..: .: :.:: ..  :  :  .:.:::::.::::::.... .::.. :.::.:.:.:..:::.. .... . ::   ::     :   . :.  .. :   ::.: .:.:: ...    .:  .: ...:.::.:.::....:: :.. .:.:..:..:  :..:: . :..  .  ..: .:   :.. .: :. ::  ..');
	}
    }
    last if( $count++ > @valid );
} 

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
ok($result->get_statistic('X1'), 6);
ok($result->get_statistic('X2'), 15);
ok($result->get_statistic('S1'), 12);
ok($result->get_statistic('S2'), 17);

ok($result->get_statistic('dbentries'), 1083200);

@valid = ( [ 'gb|AY052359.1|', 2826, 'AY052359.1', '3e-18', 96],
	   [ 'gb|AC002329.2|AC002329', 76170, 'AC002329', '3e-18', 96],
	   [ 'gb|AF132318.1|AF132318', 5383, 'AF132318', '0.040', 42]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 60);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->start, 154);
	    ok($hsp->hit->end, 212);
	    ok($hsp->hit->strand, 1);
	    ok($hsp->length('hsp'), 60);	    
	    ok($hsp->evalue == '3e-18');
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
	}
    }
    last if( $count++ > @valid );
} 

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
ok($result->query_accession, 'BACDNAE');
ok($result->query_length, 2001);
ok($result->get_parameter('matrix'), 'blosum62');

ok($result->get_statistic('lambda'), 0.318);
ok($result->get_statistic('kappa'), 0.135);
ok($result->get_statistic('entropy'),0.401 );

ok($result->get_statistic('dbentries'), 4289);

@valid = ( [ 'gi|1789447|gb|AAC76102.1|', 581, 'AAC76102.1', '1.1e-74', 671]);
$count = 0;

while( my $hit = $result->next_hit ) {
    my $d = shift @valid;
    ok($hit->name, shift @$d);
    ok($hit->length, shift @$d);
    ok($hit->accession, shift @$d);
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 21);
	    ok($hsp->query->end, 1265);
	    ok($hsp->query->strand, 1);
	    ok($hsp->hit->start, 1);
	    ok($hsp->hit->end, 413);
	    ok($hsp->hit->strand, 0);
	    ok($hsp->length('hsp'), 421);	    
	    ok($hsp->evalue == '1.1e-74');
	    ok($hsp->pvalue == '1.1e-74');
	    ok($hsp->score,671);
	    ok($hsp->bits,265.8);
	    ok(sprintf("%.2f",$hsp->percent_identity), 35.87);
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.1213);	    
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.3656);
	    ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.5297);
	    ok($hsp->query->frame(), 2);
	    ok($hsp->hit->frame(), 0);
	    ok($hsp->gaps('query'), 6);
	    ok($hsp->gaps('hit'), 8);
	    ok($hsp->gaps, 14);	    
	    ok($hsp->query_string, 'MGNRIPDEIVDQVQKSADIVEVIGDYVQLKKQGRNYFGLCPFHGESTPSFSVSPDKQIFHCFGCGAGGNVFSFLRQMEGYSFAESVSHLADKYQIDFPDDITVHSGARP---ESSGEQKMAEAHELLKKFYHHLLINTKEGQEALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGYFDRFRNRVMFPIHDHHGAVVAFSGRALGSQQPKYMNSPETPLFHKSKLLYNFYKARLHIRKQERAVLFEGFADVYTAVSSDVKESIATMGTSLTDDHVKILRRNVEEIILCYDSDKAGYEATLKASELL---QKKGCKVRVAMIPDGLDPDDYIKKFGGEKFKNDIIDASVTVMAFKMQYFRKGKNLSDEGDRLAYIKDVLKEISTLSGSLEQEVYVKQ');
	    ok($hsp->hit_string, 'MAGRIPRVFINDLLARTDIVDLIDARVKLKKQGKNFHACCPFHNEKTPSFTVNGEKQFYHCFGCGAHGNAIDFLMNYDKLEFVETVEELAAMHNLEVPFE----AGSGPSQIERHQRQTLYQLMDGLNTFYQQSL-QQPVATSARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY-DRFRERVMFPIRDKRGRVIGFGGRVLGNDTPKYLNSPETDIFHKGRQLYGLYEAQQDNAEPNRLLVVEGYMDVVALAQYGINYAVASLGTSTTADHIQLLFRATNNVICCYDGDRAGRDAAWRALETALPYMTDGRQLRFMFLPDGEDPDTLVRKEGKEAFEARM-EQAMPLSAFLFNSLMPQVDLSTPDGRARLSTLALPLISQVPGETLR-IYLRQ');
	    ok($hsp->homology_string, 'M  RIP   ++ +    DIV++I   V+LKKQG+N+   CPFH E TPSF+V+ +KQ +HCFGCGA GN   FL   +   F E+V  LA  + ++ P +    +G+ P   E    Q + +  + L  FY   L        A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y DRFR RVMFPI D  G V+ F GR LG+  PKY+NSPET +FHK + LY  Y+A+    +  R ++ EG+ DV       +  ++A++GTS T DH+++L R    +I CYD D+AG +A  +A E        G ++R   +PDG DPD  ++K G E F+  + + ++ + AF         +LS    R       L  IS + G   + +Y++Q');
	}
    }
    last if( $count++ > @valid );
} 

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
    ok($hit->significance, shift @$d );
    ok($hit->raw_score, shift @$d );

    if( $count == 0 ) {
	while( my $hsp = $hit->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 415);
	    ok($hsp->query->strand, 0);
	    ok($hsp->hit->start, 4778);
	    ok($hsp->hit->end, 6016);
	    ok($hsp->hit->strand, 1);
	    ok($hsp->length('hsp'), 421);	    
	    ok($hsp->evalue == '1.4e-73');
	    ok($hsp->pvalue == '1.4e-73');
	    ok($hsp->score,671);
	    ok($hsp->bits,265.8);
	    ok(sprintf("%.2f",$hsp->percent_identity), 35.87);
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.1219);	    
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
	}
    }
    last if( $count++ > @valid );
}

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
ok($result->query_accession, 'BACDNAE');
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
    ok($hit->significance, shift @$d );
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
		ok($hsp->evalue == '6.4e-70');
		ok($hsp->pvalue == '6.4e-70');
		ok($hsp->score,85);
		ok($hsp->bits,41.8);
		ok(sprintf("%.2f",$hsp->percent_identity), '32.20');
		ok(sprintf("%.4f",$hsp->frac_identical('hit')), 0.1073);
		ok(sprintf("%.4f",$hsp->frac_identical('query')), 0.1073);
		ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), 0.4746);
		ok($hsp->query->frame(), 2);
		ok($hsp->hit->frame(), 1);
		ok($hsp->gaps('query'), 0);
		ok($hsp->gaps('hit'), 0);
		ok($hsp->gaps, 0);	    
		ok($hsp->query_string, 'ALDYLLSRGFTKELINEFQIGYALDSWDFITKFLVKRGFSEAQMEKAGLLIRREDGSGY');
	    ok($hsp->hit_string, 'ARQYLEKRGLSHEVIARFAIGFAPPGWDNVLKRFGGNPENRQSLIDAGMLVTNDQGRSY');
	    ok($hsp->homology_string, 'A  YL  RG + E+I  F IG+A   WD + K       +   +  AG+L+  + G  Y');
	    }
	} 
    } elsif( $count == 1 ) {
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
	    ok(sprintf("%.4f",$hsp->frac_identical('hit')), '0.1250');
	    ok(sprintf("%.4f",$hsp->frac_identical('query')), '0.1250');
	    ok(sprintf("%.4f",$hsp->frac_conserved('hsp')), '0.4750');
	    ok($hsp->query->frame(), 2);
	    ok($hsp->hit->frame(), 2);
	    ok($hsp->gaps('query'), 0);
	    ok($hsp->gaps('hit'), 0);
	    ok($hsp->gaps, 0);
	    ok($hsp->query_string, 'WLPRALPEKATTAP**SWIGNMTRFLKRSKYPLPSSRLIR');
	    ok($hsp->hit_string, 'WLSRTTVGSSTVSPRTFWITRMKVKLSSSKVTLPSTKSTR');
	    ok($hsp->homology_string, 'WL R     +T +P   WI  M   L  SK  LPS++  R');
	    last;
	}       
    }
    last if( $count++ > @valid );
}

# Do a multiblast report test
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','multi_blast.bls'));

my @expected = qw(CATH_RAT CATL_HUMAN CATL_RAT PAPA_CARPA);
while( my $result = $searchio->next_result ) {
    ok($result->query_name, shift @expected, "Multiblast query test");
}

# TODO: Flesh this test out!
$searchio = new Bio::SearchIO ('-format' => 'psiblast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

$result = $searchio->next_result;

ok($result);
$hit = $result->next_hit;
ok($hit->accession, 'AE000208');
$hsp = $hit->next_hsp;
ok($hsp->get_aln->isa('Bio::Align::AlignI'));
my $writer = Bio::SearchIO::Writer::HitTableWriter->new( 
                                  -columns => [qw(
                                                  query_name
                                                  query_length
                                                  hit_name
                                                  hit_length
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
ok(-e "searchio.html");

END { 
    unlink 'searchio.out';
    unlink 'searchio.html';

}
