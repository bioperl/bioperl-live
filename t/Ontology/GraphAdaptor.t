# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 28,
			   -requires_module => 'Graph');
	
    use_ok('Bio::Ontology::SimpleGOEngine::GraphAdaptor');
}

my $g=Bio::Ontology::SimpleGOEngine::GraphAdaptor->new();
my $graph_version=( defined($Graph::VERSION) && $Graph::VERSION >= 0.5 )  ? 'new' : 'old';
my $adaptor_class=$graph_version eq 'new' ? 
  'Bio::Ontology::SimpleGOEngine::GraphAdaptor' : 'Bio::Ontology::SimpleGOEngine::GraphAdaptor02';
is (ref $g, $adaptor_class);

$g->add_vertex('vertex0');
ok($g->has_vertex('vertex0'));
ok(!$g->has_vertex('vertex1'));
my @v=$g->vertices;
is (@v, 1);
is ($v[0],'vertex0');

$g->add_edge('vertex0','vertex1');
ok($g->has_edge('vertex0','vertex1'));
ok(!$g->has_edge('vertex0','vertex'));
my @e=$g->edges;
is (@e, 1);
is ($e[0]->[0],'vertex0');
is ($e[0]->[1],'vertex1');

@e=$g->edges_at('vertex0');
is (@e, 1);
is ($e[0]->[0], 'vertex0');
is ($e[0]->[1], 'vertex1');

@v=$g->predecessors('vertex1');
is (@v, 1);
is ($v[0],'vertex0');

@v=$g->successors('vertex0');
is(@v, 1);
is ($v[0],'vertex1');

@v=$g->source_vertices;
is(@v,1);
is($v[0],'vertex0');

@v=$g->sink_vertices;
is(@v, 1);
is($v[0],'vertex1');

$g->set_vertex_attribute('vertex0','ATTR0','vertex0_ATTR0');
$g->set_vertex_attribute('vertex0','ATTR1','vertex0_ATTR1');
$g->set_vertex_attribute('vertex1','ATTR0','vertex1_ATTR0');
$g->set_vertex_attribute('vertex1','ATTR1','vertex1_ATTR1');
is ($g->get_vertex_attribute('vertex0','ATTR0'),'vertex0_ATTR0');
is ($g->get_vertex_attribute('vertex0','ATTR1'),'vertex0_ATTR1');
is ($g->get_vertex_attribute('vertex1','ATTR0'),'vertex1_ATTR0');
is ($g->get_vertex_attribute('vertex1','ATTR1'),'vertex1_ATTR1');

$g->set_edge_attribute('vertex0','vertex1','ATTR0','vertex0_vertex1_ATTR0');
$g->set_edge_attribute('vertex0','vertex1','ATTR1','vertex0_vertex1_ATTR1');
is ($g->get_edge_attribute('vertex0','vertex1','ATTR0'),'vertex0_vertex1_ATTR0');
is ($g->get_edge_attribute('vertex0','vertex1','ATTR1'),'vertex0_vertex1_ATTR1');   
