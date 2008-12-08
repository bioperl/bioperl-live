# $Id$
# -*-Perl-*-
use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;

  test_begin(-tests => 30,
	     -requires_modules => [qw(Graph::Directed
				      Bio::Tree::Node
				      Bio::TreeIO
				      IO::String
				      Array::Compare
				      Algorithm::Munkres
				      Bio::PhyloNetwork::muVector)]);

  use_ok('Bio::PhyloNetwork');
}

my $net = Bio::PhyloNetwork->new(-eNewick=>'((#H1,(#H2,l2)T),#H2)r; ((#H3,l1))#H1; ((#H3,(l3,#H1)))#H2; (l4)#H3;');

isa_ok($net, 'Bio::PhyloNetwork');

ok($net->is_leaf('l1'),'is_leaf');
ok($net->is_tree_node('T'),'is_tree_node');
ok($net->is_hybrid_node('#1'),'is_hybrid_node');
ok($net->is_root('r'),'is_root');
ok($net->is_tree_child(),'is_tree_child');

is(scalar $net->nodes(),13,'nodes');
is(scalar $net->leaves(),4,'leaves');
is(scalar $net->internal_nodes(),9,'internal_nodes');
is(scalar $net->edges(),15,'edges');
is(scalar $net->tree_edges(),9,'tree_edges');
is(scalar $net->hybrid_edges(),6,'hybrid_edges');

my %mudata = $net->mudata();
ok(($mudata{'r'} cmp 
    Bio::PhyloNetwork::muVector->new([3,1,2,5]))==0,'mudata');

my $net2 = Bio::PhyloNetwork->new(-mudata=>\%mudata,
				  -leaves=>[qw(l1 l2 l3 l4)]);
is($net->mu_distance($net2),0,'build from mudata');

my %heights=$net->heights;
is($heights{'r'},9,'heights');

my $net3 = Bio::PhyloNetwork->new(-eNewick=>'r=(l2,((l1,(#H1,l4)),#H1)); #H1=(l3);');
is($net->mu_distance($net3),12,'mu_distance');

my $net4 = Bio::PhyloNetwork->new(-eNewick=>'((1,#H1)t1,(#H1,3))r; (2)#H1');
ok($net4->is_time_consistent(),'time_consistent');

my %time=$net4->temporal_representation();
ok($time{'t1'} == $time{'#1'},'temporal_representation');

my $net5 = Bio::PhyloNetwork->new(-eNewick=>'((1,#H1)x,(#H1,3)); (2)#H1;');
my %tripartitions=$net5->tripartitions();
ok($tripartitions{'x'}  eq 'ABC','tripartitions');
ok($tripartitions{'#1'} eq 'CAC','tripartitions');

my $net6=Bio::PhyloNetwork->new(-eNewick=>'(1,(2,3))');
ok(abs($net5->tripartition_error($net6)-0.2678)<0.0001,'tripartition_error');

my ($w,$alignR,$weightsR)=$net->optimal_alignment($net3);
my %align=%{$alignR};
my %weights=%{$weightsR};

is($w,8,'weight (manhattan)');
ok($align{'T'} eq 'r','optimal_alignment (manhattan)');
my ($w2,$alignR2,$weightR2)= $net->optimal_alignment($net3,
						     -metric=>'Hamming');
is($w2,6.125,'weight (hamming)');

my %align2=%$alignR2;
ok($align2{'T'} eq 'r','optimal_alignment (hamming)');
my $enw=$net->eNewick();
my $enwf=$net->eNewick_full();


$net5=Bio::PhyloNetwork->new(-eNewick=>$enw);
ok($net->mu_distance($net5) == 0,'eNewick');

$net6=Bio::PhyloNetwork->new(-eNewick=>$enwf);
ok($net->mu_distance($net6) == 0,'eNewick_full');

my @trees=$net->explode();
is(scalar @trees,8,'explode');
isa_ok($trees[0],'Bio::Tree::Tree');
