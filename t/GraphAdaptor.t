# -*-Perl-*-
## Bioperl Test Harness Script for Modules
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 28;
}

SKIP: {
    eval { require 'Graph.pm' };
    if( $@ ) {
	    skip("Graph.pm doesn't seem to be installed on this system -- the GO Parser needs it...", 28);
    }
    use_ok('Bio::Ontology::SimpleGOEngine::GraphAdaptor');
    
    my $g=new Bio::Ontology::SimpleGOEngine::GraphAdaptor;
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
}