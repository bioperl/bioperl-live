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
    $NTESTS = 2;
    $LASTXMLTEST = 2;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 

    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	$SKIPXML = 1;
	print STDERR "XML::Parser::PerlSAX not loaded. This means ClusterIO::dbsnp test cannot be executed. Skipping\n";
	foreach ( 1..$LASTXMLTEST ) {
	    skip('No XML::Parser::PerlSAX loaded',1);
	}
    }
}

if( $error == 1 ) {
    exit(0);
}

use Bio::ClusterIO;
use Bio::Root::IO;

ok(1);
my ($clusterio, $result,$hit,$hsp);
if( ! $SKIPXML ) {
	$clusterio = new Bio::ClusterIO ('-tempfile' => 0,
					'-format' => 'dbsnp',
					'-file'   => Bio::Root::IO->catfile('t','data','LittleChrY.dbsnp.xml'));
    
	$result = $clusterio->next_cluster;
	ok($result);    
	ok($result->observed eq 'C/T');
	ok($result->type eq 'notwithdrawn');
	ok($result->seq_5);
	ok($result->seq_3);
	my @ss = $result->each_subsnp;
	ok(scalar @ss == 2);
	ok($ss[0]->handle eq 'OEFNER');
	ok($ss[1]->handle eq 'ALLENDAY');

#    ok($result->database_name, '/data_2/jason/db/cdd/cdd/Pfam');
#    ok($result->query_name,'gi|1786182|gb|AAC73112.1| (AE000111) thr operon leader peptide [Escherichia coli]');
#    ok($result->query_length, 21);
#    ok($result->algorithm, 'BLASTP');
#    ok($result->algorithm_version, 'blastp 2.1.3 [Apr-1-2001]');

#    ok($result->available_parameters, 8);
#    ok($result->get_parameter('gapext'), 1);
#    ok($result->available_statistics, 5);
#    ok($result->get_statistic('lambda'), 0.267);

# this result actually has a hit
#    $result = $searchio->next_result;
#    $hit = $result->next_hit;
#    ok($hit->name, 'gnl|Pfam|pfam00742');
#    ok($hit->description(), 'HomoS_dh, HomoS dehydrogenase');
#    ok($hit->accession, 'pfam00742');
#    ok($hit->length, 310);

#    $hsp = $hit->next_hsp;
#    ok($hsp->pvalue, undef);
#    ok($hsp->evalue, 1.46134e-90);
#    ok($hsp->score, 838);
#    ok($hsp->bits,327.405);
#    ok($hsp->query->start, 498);
#    ok($hsp->query->end,815);
#    ok($hsp->hit->start, 3);
#    ok($hsp->hit->end, 310);
#    ok($hsp->query->frame,0);
#    ok($hsp->hit->frame,0);
#    ok(sprintf("%.2f", $hsp->percent_identity), 37.73);
#    ok(sprintf("%.4f", $hsp->frac_identical('hit')), 0.3994);
#    ok(sprintf("%.4f", $hsp->frac_identical('query')), 0.3868);
#    ok(sprintf("%.4f",$hsp->query->frac_identical), 0.3868);

#    while( $result = $searchio->next_result ) { ok($result); }

}
