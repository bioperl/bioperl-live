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
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use vars qw($NTESTS);
    $NTESTS = 13;

    use Test::More;
    plan tests => $NTESTS; 
	use_ok('Bio::ClusterIO');
	use_ok('Bio::Root::IO');
	use_ok('Bio::Cluster::ClusterFactory');
}

SKIP: {
    eval { require XML::Parser::PerlSAX; };
    if( $@ ) {
		skip("XML::Parser::PerlSAX not loaded.  This means ".
			 "ClusterIO::dbsnp test cannot be executed.",8);
	}
	my ($clusterio, $result,$hit,$hsp);
	$clusterio = new Bio::ClusterIO ('-tempfile' => 0,
					'-format' => 'dbsnp',
					'-file'   => Bio::Root::IO->catfile('t','data','LittleChrY.dbsnp.xml'));
    
	$result = $clusterio->next_cluster;
	ok($result);
	is($result->observed, 'C/T');
	is($result->type, 'notwithdrawn');
	ok($result->seq_5);
	ok($result->seq_3);
	my @ss = $result->each_subsnp;
	is scalar @ss,  5;
	is($ss[0]->handle, 'CGAP-GAI');
	is($ss[1]->handle, 'LEE');
	
	# don't know if these were ever meant to work... cjf 3/7/07
	#is($result->heterozygous, 0.208738461136818);
	#is($result->heterozygous_SE, 0.0260274689436777);
}

###################################
# ClusterFactory tests            #
###################################

my $fact = Bio::Cluster::ClusterFactory->new();
# auto-recognize implementation class
my $clu = $fact->create_object(-display_id => 'Hs.2');
isa_ok($clu, "Bio::Cluster::UniGeneI");
$clu = $fact->create_object(-namespace => "UNIGENE");
isa_ok($clu, "Bio::Cluster::UniGeneI");
