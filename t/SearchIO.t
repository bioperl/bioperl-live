# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

use strict;
BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }

    use Test;
    plan tests => 47; 

    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	print STDERR "XML::Parser::PerlSAX not loaded. This means SearchIO::blastxml test cannot be executed. Skipping\n";
	foreach ( 1..43 ) {
	    skip(1,1);
	}
       $error = 1;
	
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
