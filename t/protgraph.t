# This is -*-Perl-*- code#
# Bioperl Test Harness Script for Modules#
# $Id: protgraph.t,v 1.1 2004/03/13 23:45:32 radams Exp
use vars qw($NUMTESTS $DEBUG $ERROR $XML_ERROR);
use strict;
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test;};
    $ERROR = 0;
    if ( $@ ) {
        use lib 't';
    }
    use Test;
    $NUMTESTS  = 35;
    plan tests => $NUMTESTS;
    eval {	require Class::AutoClass;	
         	require Clone;
            };
    if ( $@ ) {
        warn("Class::AutoClass or Clone not installed. " .
             " This means that the module is not usable. Skipping tests");
        $ERROR = 1;
    }
	
    eval {
        require XML::Twig;
    };
    if ($@) {
        warn "XML::Twig needed for XML format parsing, skipping these tests";
        $XML_ERROR = 1;
    }
}

END {
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
ok my $io = Bio::Graph::IO->new(-format => 'dip',
                                -file   => Bio::Root::IO->catfile("t","data","tab1part.mif"),
                                -threshold => 0.6);

ok  $gr = $io->next_network();
#ok my @conns = $gr->components();
#print STDERR scalar @conns, "s\n";
 my @nodes = $gr->articulation_points();
 for my $node(@nodes){
		my $n   = $gr->_nodes($node);
		my $ac  = $n->annotation();
		my @dbs = $ac->get_Annotations('dblink');
		my $id  = (map{$_->primary_id}grep{$_->database eq 'DIP'} @dbs)[0];
		$ids{$id}++;
	}

for (sort keys %ids) {
	print "$_ : $ids{$_}\n";
}		
ok my $node   = $gr->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

##test write. to filehandle...

ok my $out =  Bio::Graph::IO->new(-format => 'dip',
                                  -file   =>">". Bio::Root::IO->catfile("t","data","out.mif"));
ok $out->write_network($gr);

## can we round trip, is out format same as original format?
ok my $io2 = Bio::Graph::IO->new(-format    => 'dip',
                                 -file     => Bio::Root::IO->catfile("t","data","out.mif"));
ok	my $g2     = $io2->next_network();
ok  $node      = $g2->nodes_by_id('A64696');
ok $node->accession_number, 'A64696';

############## now lets test the XML format.....
if (!$XML_ERROR){


}


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
my $i = 0;
for my $n ($g2->nodes) {
 if ( ref($n) !~ /Bio/) {
		print "-$i- - heere";
		}
$i++;
};

## count all edges
my $count = 0;
ok scalar keys %{$g2->_edges}, 71;

my @n = $g2->neighbors($g2->nodes_by_id('3075N'));
ok scalar @n, 13;

ok $g2->remove_nodes($g2->nodes_by_id('3075N'));

## should be 13  less interactions in graph.  
ok scalar keys %{$g2->_edges}, 58;

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

