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
    $NTESTS = 10;
    $LASTXMLTEST = 10;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 

    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
	$SKIPXML = 1;
	print STDERR "XML::Parser::PerlSAX not loaded. This means ClusterIO::dbsnp test cannot be executed. Skipping\n";
	foreach ( $Test::ntest..$LASTXMLTEST ) {
	    skip('No XML::Parser::PerlSAX loaded',1);
	}
    }
}

if( $error == 1 ) {
    exit(0);
}

use Bio::ClusterIO;
use Bio::Root::IO;

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
	ok($result->heterozygous == 0.208738461136818);
	ok($result->heterozygous_SE == 0.0260274689436777);
}
