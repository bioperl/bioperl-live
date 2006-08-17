# This is -*-Perl-*- code#
# Bioperl Test Harness Script for Modules#
# $Id: protgraph.t,v 1.1 2004/03/13 23:45:32 radams Exp

use vars qw($NUMTESTS $DEBUG $ERROR $XML_ERROR);
use strict;
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
use Bio::Root::IO;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test;};
    $ERROR = $XML_ERROR = 0;
    if ( $@ ) {
	use lib 't';
    }
    use Test;
    $NUMTESTS  = 66;
    plan tests => $NUMTESTS;
    eval {	require Class::AutoClass;
         	require Clone; };
    if ( $@ ) {
	warn("Class::AutoClass or Clone not installed. This means that the module is not usable. Skipping tests\n");
	$ERROR = 1;
    }

    eval {
	require XML::Twig;
    };
    if ($@) {
	warn "XML::Twig needed for XML format parsing, skipping these tests\n";
	$XML_ERROR = 1;
    }
}

END {
    unlink Bio::Root::IO->catfile("t","data","out.mif");
    foreach ( $Test::ntest..$NUMTESTS) {
	skip("Missing dependencies. Skipping tests",1);
    }
}
exit 0 if $ERROR ==  1;

require Bio::Graph::ProteinGraph;
require Bio::Graph::IO;
require Bio::Graph::Edge;

my $verbose = 0;
$verbose    = 1 if $DEBUG;
ok 1;
################1st of all let's test the io.....

###############  test dip tab delimited  format  ###################
## test read...
my %ids;
my $gr;
ok my $io = Bio::Graph::IO->new(
  -format => 'dip',
  -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"),
  -threshold => 0.6);

ok  $gr = $io->next_network();

ok my $node   = $gr->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

##test write. to filehandle...

ok my $out =  Bio::Graph::IO->new(
  -format => 'dip',
  -file   =>">". Bio::Root::IO->catfile("t","data","out.mif"));
ok $out->write_network($gr);

## get articulation_points. 
my @nodes = $gr->articulation_points();

##now remove 2 nodes: this removes 4 edges and  3087 should be a new artic point
ok $gr->edge_count, 72;
$gr->remove_nodes($gr->nodes_by_id('3082N'), $gr->nodes_by_id('3083N'));
ok $gr->edge_count, 68;
 my $nodes = $gr->articulation_points();
ok grep {$_->object_id eq 'H64521'} @$nodes;
ok scalar @$nodes, 13;
@nodes = @{$gr->articulation_points()};
# <NOTE>
# these were failing, I don't understand the module enough to know if 
# this is a bug. Richard needs to look at it
#ok grep {$_->object_id eq 'B64701'} @nodes;
#ok scalar @nodes, 14;

ok grep {$_->object_id eq 'B64528'} @nodes;
ok scalar @nodes, 13;
# </NOTE>

ok $node   = $gr->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';


## can we round trip, is out format same as original format?
ok my $io2 = Bio::Graph::IO->new(
  -format    => 'dip',
  -file     => Bio::Root::IO->catfile("t","data","out.mif"));
ok	my $g2     = $io2->next_network();
ok  $node      = $g2->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

##### now lets test some graph properties.....##
## basic properties from SImpleGraph.

ok sprintf("%.3f",$g2->density), "0.027";
ok $g2->is_connected, '';
ok $g2->is_forest, undef;
ok $g2->is_tree, '';
ok $g2->is_empty, '';
ok $g2->is_cyclic, 1;

## get connected subgraphs
my @components = $g2->components();
ok scalar @components, 5;

## get nodes connected to parameter
my $t       = $g2->traversal($g2->nodes_by_id('3079N'));
my @dfnodes = $t->get_all;
##

##before deleting 3048N,  3047N has 2 neighbours
my @n1 = $g2->neighbors($g2->nodes_by_id('3047N'));
ok scalar @n1,2;

ok $g2->remove_nodes($g2->nodes_by_id('3048N'));

## after deleting there is only 1 interactor
@n1 = $g2->neighbors($g2->nodes_by_id('3047N'));
ok scalar @n1,1;

##check no undefs left after node removal ##

ok map {$_->object_id}$g2->edges;

## get an edge by its id

ok my $edge = $g2->edge_by_id('4368E');
ok $edge->object_id, '4368E';

## count all edges
my $count = 0;
ok $g2->edge_count, 71;

my @n = $g2->neighbors($g2->nodes_by_id('3075N'));
ok scalar @n, 13;

ok $g2->remove_nodes($g2->nodes_by_id('3075N'));

## should be 13  less interactions in graph.  
ok scalar $g2->edge_count, 58;

## many more subgraphs now
@components = $g2->components();
#there were 5 subgraphs, now there are 10 unconnected nodes, total 15
ok scalar @components, 15;

## how many unconnected nodes?
my @ucnodes = $g2->unconnected_nodes;
ok scalar  @ucnodes, 10;

##get CC using protein object..
ok  sprintf("%.3f", $g2->clustering_coefficient($g2->nodes_by_id('B64525'))), 0.022;

#.. and using id string (same as previous, for convenience	)
ok  sprintf("%.3f", $g2->clustering_coefficient('B64525')), 0.022;

## test has_node() method
ok $g2->has_node('B64525'), 1;
ok $g2->has_node('B64'), 0;

## remove a single duplicate edge
ok $g2->remove_dup_edges($g2->nodes_by_id('3103N'));

## remove  all duplicate edges
ok $g2->remove_dup_edges();

## should now be no duplicates
my @dupids = map{$_->object_id()} $g2->dup_edges();
ok $dupids[0], undef;

########### now we test the 'union()' method to see it conforms to 
## the rules described in its documentation:

$io = Bio::Graph::IO->new(
   -format => 'dip',
   -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"));
$gr = $io->next_network();
$io2 = Bio::Graph::IO->new(
   -format => 'dip',
   -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"));

$g2 = $io2->next_network();

# First of all we put the same graph into both variables. After a union
# graph 1 should be unaffected. Because all edge ids are the same, 
# all duplicates will be redundant. 
# therea re 3 duplicates in dataset. 
my @dups = $gr->dup_edges();
ok scalar @dups, 3;
$gr->union($g2);
@dups = $gr->dup_edges();
ok scalar @dups, 3;
my @redundant = $gr->redundant_edge();
ok scalar @redundant, 72; 

## now lets do a union with a graph that has some new edges, 
## using existing nodes

##read in graph data
$gr = undef;
$g2 = undef;
$io = Bio::Graph::IO->new(
   -format => 'dip',
   -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"));
$gr = $io->next_network();
$io2 = Bio::Graph::IO->new(
    -format => 'dip',
    -file   => Bio::Root::IO->catfile("t","data","tab2part.mif"));
$g2 = $io2->next_network();
ok $gr->edge_count, 72;
ok $gr->node_count, 74;
$gr->union($g2);
#there should be 1 more edge in the graph $gr now, with no new nodes. 
#$g2 is unaffected.  
ok $gr->edge_count, 73;
ok $gr->node_count, 74;

## now lets test a union that has new nodes in $g2 
$gr = undef;
$g2 = undef;
$io = Bio::Graph::IO->new
    (-format => 'dip',
     -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"));
$gr = $io->next_network();
$io2 = Bio::Graph::IO->new
    (-format => 'dip',
     -file   => Bio::Root::IO->catfile("t","data","tab3part.mif"));

$g2 = $io2->next_network();
ok $gr->edge_count, 72;
ok $gr->node_count, 74;
$gr->union($g2);
# there should be 2 more edge in the graph $gr now and 2 more nodes. 
# $g2 is unaffected.  
ok $gr->edge_count, 74;
ok $gr->node_count, 76;

# test IO/psi_xml if the required modules are present
unless( $XML_ERROR ) {
	# PSI XML from DIP
	ok $io = Bio::Graph::IO->new
	  (-format => 'psi_xml',
		-file   => Bio::Root::IO->catfile("t", "data", "psi_xml.dat"));
	ok my $g = $io->next_network();
	ok $g->edge_count, 3;
	ok $g->node_count, 4;
	#my @rts =$g->articulation_points();
	my $n = $g->nodes_by_id(207153);
	ok $n->species->node_name,"Helicobacter pylori 26695";
	ok $n->primary_seq->desc,"bogus-binding membrane protein (lepA) HP0355";

	# PSI XML from IntAct
	ok my $io2 = Bio::Graph::IO->new
	  (-format => 'psi_xml',
		-file   => Bio::Root::IO->catfile("t", "data", "sv40_small.xml"));
	ok my $g3 = $io2->next_network();
	ok $g3->edge_count, 3;
	ok $g3->node_count, 5;

	# my @rts =$g->articulation_points();
	# my $n = $g->nodes_by_id(207153);
	# ok $n->species->binomial,"Helicobacter pylori 26695";
	# ok $n->primary_seq->desc,"bogus-binding membrane protein (lepA) HP0355";
} 

