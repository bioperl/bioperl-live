# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use vars qw($NTESTS);
    $NTESTS = 245;
    use Test;
    plan tests => $NTESTS; 

    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	print STDERR "XML::Parser::PerlSAX not loaded. This means SearchIO::blastxml test cannot be executed. Skipping\n";
	foreach ( 1..$NTESTS ) {
	    skip(1,1);
	}
       $error = 1;	
    } else {
	$error = 0;
    }
}

if( $error == 1 ) {
    exit(0);
}

use Bio::SearchIO;
use Bio::Root::IO;

ok(1);

# test with RPSBLAST data first 
my $searchio = new Bio::SearchIO ('-format' => 'blastxml',
     '-file'   => Bio::Root::IO->catfile('t','data','ecoli_domains.rps.xml'));

my $report = $searchio->next_report;
ok($report);
ok($report->database_name, '/data_2/jason/db/cdd/cdd/Pfam');
ok($report->query_name,'gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide [Escherichia coli]');
ok($report->query_size, 21);
ok($report->program_name, 'blastp');
ok($report->program_version, 'blastp 2.1.3 [Apr-1-2001]');

ok($report->available_parameters, 8);
ok($report->get_parameter('gapext'), 1);
ok($report->available_statistics, 7);
ok($report->get_statistic('lambda'), 0.267);

# this report actually has a hit
$report = $searchio->next_report;
my $subject = $report->next_subject;
ok($subject->name, 'gnl|Pfam|pfam00742');
ok($subject->description, 'HomoS_dh, HomoS dehydrogenase');
ok($subject->accession, 'pfam00742');
ok($subject->length, 310);

my $hsp = $subject->next_hsp;
ok($hsp->P, 1.46134e-90);
ok($hsp->evalue, 1.46134e-90);
ok($hsp->score, 838);
ok($hsp->bits,327.405);
ok($hsp->query->start, 498);
ok($hsp->query->end,815);
ok($hsp->subject->start, 3);
ok($hsp->subject->end, 310);
ok($hsp->frame,0);
ok($hsp->hframe,0);
ok($hsp->percent_identity, );
ok($hsp->subject->frac_identical, 123);
ok($hsp->query->frac_identical, 123);

ok($hsp->positive,171);
ok($hsp->gaps, 26);
ok($hsp->hsp_length, 326);
ok($hsp->query_seq, 'LRVCGVANSKALLTNVHGLNLENWQEELAQAKEPF-NLGRLIRLVKEYHLLN----PVIVDCTSSQAVAD-QYADFLREGFHVVTPNKKANTSSMDYYHQLRYAAEKSRRKFLYDTNVGAGLPVIENLQNLLNAGDELMKFSGILSGSLSYIFGKLDE-GMSFSEATTLAREMGYTEPDPRDDLSGMDVARKLLILARET-GRELELADIEIEPVLPAEFNAEGDVAAFMANLSQLDDLFAARVAKARDEGKVLRYVGNIDEDGVCRVKIAEVDGNDPLFKVKNGENALAFYSHYYQPLPLVLRGYGAGNDVTAAGVFADLLRTLS');
ok($hsp->subject_seq, 'GVVTGITDSREMLLSRIGLPLEIWKVALRDLEKPRKDLGKLDLTDDAFAVVDDPDIDVVVELTGGIEVARELYLDALEEGKHVVTANKALNASHGDEYLAL---AEKSGVDVLYEAAVAGGIPIIKTLRELLATGDRILKIEGIFNGTTNFILSEMDEKGLPFSDVLAEAQELGYTEADPRDDVEGIDAARKLAILARIAFGIELELDDVYVEGISPITAEDISSADEFGYTLKLLDEAMRQRVEDAESGGEVLRYPTLIPE-------------DHPLASVKGSDNAVAVEGEAYG--PLMFYGPGAGAEPTASAVVADIVRIAR');
ok($hsp->homology_seq, '  V G+ +S+ +L +  GL LE W+  L   ++P  +LG+L      + +++     V+V+ T    VA   Y D L EG HVVT NK  N S  D Y  L   AEKS    LY+  V  G+P+I+ L+ LL  GD ++K  GI +G+ ++I  ++DE G+ FS+    A+E+GYTE DPRDD+ G+D ARKL ILAR   G ELEL D+ +E + P           F   L  LD+    RV  A   G+VLRY   I E             + PL  VK  +NA+A     Y   PL+  G GAG + TA+ V AD++R   ');

# one more 
$subject = $report->next_subject;
ok($subject);

while( $report = $searchio->next_report ) { ok($report); }

     
$searchio = new Bio::SearchIO(-format => 'blastxml', -file => Bio::Root::IO->catfile('t','data','plague_yeast.bls.xml'));

$report = $searchio->next_report;

ok($report->database_name, 'yeast.aa');
ok($report->query_name, 'gi|5763811|emb|CAB53164.1| putative transposase [Yersinia pestis]');
ok($report->query_size, 340);

$subject = $report->next_subject;
ok(! $subject);


$searchio = new Bio::SearchIO ('-format' => 'blast',
				  '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.bls'));

$report = $searchio->next_report;

ok($report->database_name, 'ecoli.aa');
ok($report->database_size, 1358990);
ok($report->program_name, 'BLASTP');
ok($report->program_version, '2.1.3');
ok($report->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($report->query_size, 820);
ok($report->get_statistic('kappa')== 0.041);
ok($report->get_statistic('lambda'), 0.267);
ok($report->get_statistic('entropy') == 0.14);
ok($report->get_statistic('dblength'), 1358990);
ok($report->get_statistic('dbnum'), 4289);
ok($report->get_statistic('hsplength'), 47);
ok($report->get_statistic('effectivespace'), 894675611);
ok($report->get_parameter('matrix'), 'BLOSUM62');
ok($report->get_parameter('gapopen'), 11);
ok($report->get_parameter('gapext'), 1);

my @valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113'],
	      [ 'gb|AAC76922.1|', 810, 'AAC76922'],
	      [ 'gb|AAC76994.1|', 449, 'AAC76994']);
my $count = 0;
while( $subject = $report->next_subject ) {
    my $d = shift @valid;
    ok($subject->name, $d->[0]);
    ok($subject->length, $d->[1]);
    ok($subject->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $subject->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 820);
	    ok($hsp->subject->start, 1);
	    ok($hsp->subject->end, 820);
	    ok($hsp->hsp_length, 820);
	    ok($hsp->P == 0.0);
	    ok($hsp->evalue == 0.0);
	    ok($hsp->score, 4058);
	    ok($hsp->bits,1567);	    
	    ok($hsp->positive, 806);
	    ok($hsp->percent_identity, 98.2);
	    ok($hsp->query->frac_identical, 806);
	    ok($hsp->subject->frac_identical, 806);
	    ok($hsp->gaps, 0);	    
	}
    }
    last if( $count++ > @valid );
}

$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','ecolitst.wublastn'));

$report = $searchio->next_report;

ok($report->database_name, 'ecoli.aa');
ok($report->database_size, 1358990);
ok($report->program_name, 'BLASTP');
ok($report->program_version, '2.0MP-WashU');
ok($report->query_name, qr/gi|1786183|gb|AAC73113.1| (AE000111) aspartokinase I,\s+homoserine dehydrogenase I [Escherichia coli]/);
ok($report->query_size, 820);
ok($report->get_statistic('kappa'), 0.136);
ok($report->get_statistic('lambda'), 0.319);
ok($report->get_statistic('entropy'), 0.384);
ok($report->get_statistic('dblength'), 1358990);
ok($report->get_statistic('dbnum'), 4289);
ok($report->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AAC73113.1|', 820, 'AAC73113'],
	   [ 'gb|AAC76922.1|', 810, 'AAC76922'],
	   [ 'gb|AAC76994.1|', 449, 'AAC76994']);
$count = 0;
while( $subject = $report->next_subject ) {
    my $d = shift @valid;
    ok($subject->name, $d->[0]);
    ok($subject->length, $d->[1]);
    ok($subject->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $subject->next_hsp ) {
	    ok($hsp->query->start, 1);
	    ok($hsp->query->end, 820);
	    ok($hsp->subject->start, 1);
	    ok($hsp->subject->end, 820);
	    ok($hsp->hsp_length, 820);
	    ok($hsp->P == 0.0);
	    ok($hsp->evalue == 0.0);
	    ok($hsp->score, 4141);
	    ok($hsp->bits,1462.8);	    
	    ok($hsp->positive, 820);
	    ok($hsp->percent_identity, 100);
	    ok($hsp->query->frac_identical, 820);
	    ok($hsp->subject->frac_identical, 820);
	    ok($hsp->gaps, 0);	    
	}
    }
    last if( $count++ > @valid );
}

# test tblastx 
$searchio = new Bio::SearchIO ('-format' => 'blast',
			       '-file'   => Bio::Root::IO->catfile('t','data','HUMBETGLOA.tblastx'));

$report = $searchio->next_report;
ok($report->database_name, 'ecoli.nt');
ok($report->database_size, 4662239);
ok($report->program_name, 'TBLASTX');
ok($report->program_version, '2.1.2');
ok($report->query_name, qr/HUMBETGLOA Human haplotype C4 beta-globin gene, complete cds./);
ok($report->query_size, 3002);
ok($report->get_statistic('kappa'), 0.135);
ok($report->get_statistic('lambda'), 0.318);
ok($report->get_statistic('entropy'), 0.401);
ok($report->get_statistic('dblength'), 4662239);
ok($report->get_statistic('dbnum'), 400);
ok($report->get_parameter('matrix'), 'BLOSUM62');

@valid = ( [ 'gb|AE000479.1|AE000479', 10934, 'AE000479'],
	   [ 'gb|AE000302.1|AE000302', 10264, 'AE000302'],
	   [ 'gb|AE000277.1|AE000277', 11653, 'AE000277']);
$count = 0;

while( $subject = $report->next_subject ) {
    my $d = shift @valid;
    ok($subject->name, $d->[0]);
    ok($subject->length, $d->[1]);
    ok($subject->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $subject->next_hsp ) {
	    ok($hsp->query->start, 1057);
	    ok($hsp->query->end, 1134);
	    ok($hsp->query->strand, 1);
	    ok($hsp->subject->end, 5893);
	    ok($hsp->subject->start, 5816);
	    ok($hsp->subject->strand, -1);
	    ok($hsp->hsp_length, 26);
	    ok($hsp->P == 0.13);
	    ok($hsp->evalue == 0.13);
	    ok($hsp->score, 67);
	    ok($hsp->bits,33.6);	    
	    ok($hsp->positive, 16);
	    ok($hsp->percent_identity, 42.3);
	    ok($hsp->query->frac_identical, 11);
	    ok($hsp->subject->frac_identical, 11);
	    ok($hsp->query->frame(), 0);
	    ok($hsp->subject->frame(), 1);
	    ok($hsp->gaps, 0);	    
	    ok($hsp->query_seq, 'SAYWSIFPPLGCWWSTLGPRGSLSPL');
	    ok($hsp->subject_seq, 'AAVWALFPPVGSQWGCLASQWRTSPL');
	    ok($hsp->homology_seq, '+A W++FPP+G  W  L  +   SPL');
	}
    }
    last if( $count++ > @valid );
}

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/HUMBETGLOA.FASTA');
$report = $searchio->next_report;
ok($report->database_name, qr/dros_clones.2.5/);
ok($report->database_size, 112936249);
ok($report->program_name, 'FASTA');
ok($report->program_version, '3.3t08');
ok($report->query_name, qr/HUMBETGLOA Human haplotype C4 beta-globin gene, complete cds./);
ok($report->query_size, 3002);
ok($report->get_parameter('gapopen'), -16);
ok($report->get_parameter('gapext'), -4);
ok($report->get_parameter('ktup'), 6);

ok($report->get_statistic('lambda'), 0.0823);
ok($report->get_statistic('dblength'), 112936249);
ok($report->get_statistic('dbnum'), 657);

@valid = ( [ 'BACR21I23', 73982, 'BACR21I23'],
	   [ 'BACR40P19', 73982, 'BACR40P19'],
	   [ 'BACR30L17', 32481, 'BACR30L17']);
$count = 0;

while( my $subject = $report->next_subject ) {
    my $d = shift @valid;
    ok($subject->name, $d->[0]);
    ok($subject->length, $d->[1]);
    ok($subject->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $subject->next_hsp ) {
	    ok($hsp->query->start, 31);
	    ok($hsp->query->end, 289);
	    ok($hsp->query->strand, -1);
	    ok($hsp->subject->end, 65167);
	    ok($hsp->subject->start, 64902);
	    ok($hsp->subject->strand, 1);
	    ok($hsp->hsp_length, 267);
	    ok($hsp->P == 0.017);
	    ok($hsp->evalue == 0.017);
	    ok($hsp->score, 134.5);
	    ok($hsp->bits,44.2);
	    ok($hsp->percent_identity, 57.3);
	    ok($hsp->query->frac_identical, 153);  # not sure this is right
	    ok($hsp->subject->frac_identical, 153);  # not sure this is right
	    ok($hsp->query->frame(), 0);
	    ok($hsp->subject->frame(), 0);
	    ok($hsp->gaps, 159);	    
	    ok($hsp->query_seq, 'GATTAAAACCTTCTGGTAAGAAAAGAAAAAATATATATATATATATATGTGTATATGTACACACATACATATACATATATATGCATTCATTTGTTGTTGTTTTTCTTAATTTGCTCATGCATGCTA----ATAAATTATGTCTAAAAATAGAAT---AAATACAAATCAATGTGCTCTGTGCATTA-GTTACTTATTAGGTTTTGGGAAACAAGAGGTAAAAAACTAGAGACCTCTTAATGCAGTCAAAAATACAAATAAATAAAAAGTCACTTACAACCCAAAGTGTGACTATCAATGGGGTAATCAGTGGTGTCAAATAGGAGGT');
	    ok($hsp->subject_seq, 'GATGTCCTTGGTGGATTATGGTGTTAGGGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATACAAAATATAATATAAAATATAATATAAAATATAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATATAAAATAAAATAT-AATATAAAATATAAAATAAAATATAATATAAAATATAATATAAAATATAATATAAAATATAATATAAAATA');
	    ok($hsp->homology_seq, '                              :::::::::::::::::: : ::::: :: : : ::: ::::: ::::::::  ::  :: : :   : : : : :  ::    : :: ::   ::    : ::: :::     :::::: :::   ::::: ::  :::  :    :    : ::   :::  : ::   : :   : : :: :   :: : : :: : :       ::  : : ::: ::: ::  ::::: ::: : :  :: ::   ::: : : : ::: ::                                                               ');
	}
    }
    last if( $count++ > @valid );
} 

$searchio = new Bio::SearchIO(-format => 'fasta',
				 -file   => 't/data/cysprot1.FASTA');
$report = $searchio->next_report;
ok($report->database_name, qr/ecoli.aa/);
ok($report->database_size, 1358987);
ok($report->program_name, 'FASTA');
ok($report->program_version, '3.3t08');
ok($report->query_name, 'CYS1_DICDI');
ok($report->query_size, 343);
ok($report->get_parameter('gapopen'), -12);
ok($report->get_parameter('gapext'), -2);
ok($report->get_parameter('ktup'), 2);

ok($report->get_statistic('lambda'), 0.1456);
ok($report->get_statistic('dblength'), 1358987);
ok($report->get_statistic('dbnum'), 4289);


@valid = ( [ 'gi|1787478|gb|AAC74309.1|', 512, 'AAC74309'],
	   [ 'gi|1790635|gb|AAC77148.1|', 251, 'AAC77148'],
	   [ 'gi|1786590|gb|AAC73494.1|', 94, 'AAC73494']);
$count = 0;

while( my $subject = $report->next_subject ) {
    my $d = shift @valid;
    ok($subject->name, $d->[0]);
    ok($subject->length, $d->[1]);
    ok($subject->accession, $d->[2]);
    if( $count == 0 ) {
	while( my $hsp = $subject->next_hsp ) {
	    ok($hsp->query->start, 125);
	    ok($hsp->query->end, 305);
	    ok($hsp->query->strand, 1);
	    ok($hsp->subject->start, 2);
	    ok($hsp->subject->end, 181);
	    ok($hsp->subject->strand, 1);
	    ok($hsp->hsp_length, 188);
	    ok($hsp->P == 1.2);
	    ok($hsp->evalue == 1.2);
	    ok($hsp->score, 109.2);
	    ok($hsp->bits,29.2);
	    ok($hsp->percent_identity, 23.9);
	    ok($hsp->query->frac_identical, 45); # not sure this is right
	    ok($hsp->subject->frac_identical, 45); # not sure this is right
	    ok($hsp->query->frame(), 0);
	    ok($hsp->subject->frame(), 0);
	    ok($hsp->gaps, 49);	    
	    ok($hsp->query_seq, 'NKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTT-GNV----EGQHFISQNKLVSLSEQNLVDCDHECME-YEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGP-LAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII');
	    ok($hsp->subject_seq, 'MKIRSQVGMVLNLDKCIGCHTCSVTCKNVWTSREGVEYAWFNNVETKPGQGF-PTDWENQEKYKGGWI--RKINGKLQPRMGNRAMLLGKIFANPHLPGIDDYYEPFDFDYQNLHTAPEG----SKSQPIARPRSLITGERMAKIEKGPNWEDDLGGEFDKLAKDKNFDN-IQKAMYSQFENTFMMYLPRLCEHCLNPACVATCPSGAIYKREEDGIVLIDQDKCRGWRMCITGCPYKKIYFNWKSGKSEKCIFCYPRIEAGQPTVCSETC');
	    ok($hsp->homology_seq, '                              . :. :  : :  .: .: . :.:  ::    :: ..   :.. .   :..   : : .: :.:     .  :: :::   :  .  : : ..   :   .     .:.  :. .   .     :.. .     . ::  .:    . .:.  .:: ::   . ...:. :  . ::  .. :   .:                                                                                  ');
	}
    }
    last if( $count++ > @valid );
} 

