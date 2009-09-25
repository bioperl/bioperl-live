# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 12,
			   -requires_module => 'Time::HiRes');
    
	use_ok('Bio::ClusterIO');
	use_ok('Bio::Cluster::ClusterFactory');
}

SKIP: {
	test_skip(-tests => 8, -requires_module => 'XML::SAX');
    
	my ($clusterio, $result,$hit,$hsp);
	$clusterio = Bio::ClusterIO->new('-tempfile' => 0,
					'-format' => 'dbsnp',
					'-file'   => test_input_file('LittleChrY.dbsnp.xml'));
    
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
