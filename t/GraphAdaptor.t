# -*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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

    eval { require 'Graph.pm' };
    if( $@ ) {
	    print STDERR "\nGraph.pm doesn't seem to be installed on this system -- the GO Parser needs it...\n\n";
	    plan tests => 1;
	    ok( 1 );
	    exit( 0 );
    }

    plan tests => 18;
}

use Bio::Ontology::SimpleGOEngine::GraphAdaptor;

my $g=new Bio::Ontology::SimpleGOEngine::GraphAdaptor;
my $graph_version=( defined($Graph::VERSION) && $Graph::VERSION >= 0.5 )  ? 'new' : 'old';
my $adaptor_class=$graph_version eq 'new' ? 
  'Bio::Ontology::SimpleGOEngine::GraphAdaptor' : 'Bio::Ontology::SimpleGOEngine::GraphAdaptor02';
ok (ref $g, $adaptor_class);

$g->add_vertex('vertex0');
ok($g->has_vertex('vertex0'));
ok(!$g->has_vertex('vertex1'));
my @v=$g->vertices;
ok (@v==1 && $v[0] eq 'vertex0') ;

$g->add_edge('vertex0','vertex1');
ok($g->has_edge('vertex0','vertex1'));
ok(!$g->has_edge('vertex0','vertex'));
my @e=$g->edges;
ok (@e==1 && $e[0]->[0] eq 'vertex0' && $e[0]->[1] eq 'vertex1') ;

@e=$g->edges_at('vertex0');
ok (@e==1 && $e[0]->[0] eq 'vertex0' && $e[0]->[1] eq 'vertex1') ;

@v=$g->predecessors('vertex1');
ok (@v==1 && $v[0] eq 'vertex0');

@v=$g->successors('vertex0');
ok (@v==1 && $v[0] eq 'vertex1');

@v=$g->source_vertices;
ok (@v==1 && $v[0] eq 'vertex0');

@v=$g->sink_vertices;
ok (@v==1 && $v[0] eq 'vertex1');

$g->set_vertex_attribute('vertex0','ATTR0','vertex0_ATTR0');
$g->set_vertex_attribute('vertex0','ATTR1','vertex0_ATTR1');
$g->set_vertex_attribute('vertex1','ATTR0','vertex1_ATTR0');
$g->set_vertex_attribute('vertex1','ATTR1','vertex1_ATTR1');
ok ($g->get_vertex_attribute('vertex0','ATTR0'),'vertex0_ATTR0');
ok ($g->get_vertex_attribute('vertex0','ATTR1'),'vertex0_ATTR1');
ok ($g->get_vertex_attribute('vertex1','ATTR0'),'vertex1_ATTR0');
ok ($g->get_vertex_attribute('vertex1','ATTR1'),'vertex1_ATTR1');

$g->set_edge_attribute('vertex0','vertex1','ATTR0','vertex0_vertex1_ATTR0');
$g->set_edge_attribute('vertex0','vertex1','ATTR1','vertex0_vertex1_ATTR1');
ok ($g->get_edge_attribute('vertex0','vertex1','ATTR0'),'vertex0_vertex1_ATTR0');
ok ($g->get_edge_attribute('vertex0','vertex1','ATTR1'),'vertex0_vertex1_ATTR1');


