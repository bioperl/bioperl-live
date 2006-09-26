# $Id$
#
# BioPerl module for Bio::Graph::ProteinGraph
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::ProteinGraph - a representation of a protein interaction graph.

=head1 SYNOPSIS

  # Read in from file
  my $graphio = Bio::Graph::IO->new(-file   => 'myfile.dat',
                                    -format => 'dip');
  my $graph   = $graphio->next_network();

=head2 Using ProteinGraph

  # Remove duplicate interactions from within a dataset
  $graph->remove_dup_edges();

  # Get a node (represented by a sequence object) from the graph.
  my $seqobj = $gr->nodes_by_id('P12345');

  # Get clustering coefficient of a given node.
  my $cc = $gr->clustering_coefficient($graph->nodes_by_id('NP_023232'));
  if ($cc != -1) {  ## result is -1 if cannot be calculated
    print "CC for NP_023232 is $cc";
  }

  # Get graph density
  my $density = $gr->density();

  # Get connected subgraphs
  my @graphs = $gr->components();

  # Remove a node
  $gr->remove_nodes($gr->nodes_by_id('P12345'));

  # How many interactions are there?
  my $count = $gr->edge_count;

  # How many nodes are there?
  my $ncount = $gr->node_count();

  # Let's get interactions above a threshold confidence score.
  my $edges = $gr->edges;
  for my $edge (keys %$edges) {
	 if (defined($edges->{$edge}->weight()) &&
      $edges->{$edge}->weight() > 0.6) {
		    print $edges->{$edge}->object_id(), "\t",
             $edges->{$edge}->weight(),"\n";
	 }
  }

  # Get interactors of your favourite protein
  my $node      = $graph->nodes_by_id('NP_023232');
  my @neighbors = $graph->neighbors($node); 
  print "      NP_023232 interacts with ";
  print join " ,", map{$_->object_id()} @neighbors;
  print "\n";

  # Annotate your sequences with interaction info
  my @seqs; ## array of sequence objects
  for my $seq(@seqs) {
    if ($graph->has_node($seq->accession_number)) {
       my $node = $graph->nodes_by_id( $seq->accession_number);
       my @neighbors = $graph->neighbors($node);
       for my $n (@neighbors) {
         my $ft = Bio::SeqFeature::Generic->new(
                      -primary_tag => 'Interactor',
                      -tags        => { id => $n->accession_number }
                      );
            $seq->add_SeqFeature($ft);
        }
     }
  }

  # Get proteins with > 10 interactors
  my @nodes = $graph->nodes();
  my @hubs;
  for my $node (@nodes) {
    if ($graph->neighbor_count($node) > 10) {
       push @hubs, $node;
    }
  }
  print "the following proteins have > 10 interactors:\n";
  print join "\n", map{$_->object_id()} @hubs;

  # Merge graphs 1 and 2 and flag duplicate edges
  $g1->union($g2);
  my @duplicates = $g1->dup_edges();
  print "these interactions exist in $g1 and $g2:\n";
  print join "\n", map{$_->object_id} @duplicates;

=head2 Creating networks from your own data

If you have interaction data in your own format, e.g. 

  edgeid  node1  node2  score

  my $io = Bio::Root::IO->new(-file => 'mydata');
  my $gr = Bio::Graph::ProteinGraph->new();
  my %seen = (); # to record seen nodes
  while (my $l = $io->_readline() ) {

  # Parse out your data...
  my ($e_id, $n1, $n2, $sc) = split /\s+/, $l;

  # ...then make nodes if they don't already exist in the graph...
  my @nodes =();
    for my $n ($n1, $n2 ) {
		if (!exists($seen{$n})) {
        push @nodes,  Bio::Seq->new(-accession_number => $n);
		  $seen{$n} = $nodes[$#nodes];
      } else {
			push @nodes, $seen{$n};
	   }
    }
  }

  # ...and add a new edge to the graph
  my $edge  = Bio::Graph::Edge->new(-nodes => \@nodes,
                                    -id    => 'myid',
                                    -weight=> 1);
  $gr->add_edge($edge);

=head1 DESCRIPTION

A ProteinGraph is a representation of a protein interaction network.
It derives most of its functionality from the L<Bio::Graph::SimpleGraph>
module, but is adapted to be able to use protein identifiers to
identify the nodes.

This graph can use any objects that implement L<Bio::AnnotatableI> and 
L<Bio::IdentifiableI> interfaces.  L<Bio::Seq> (but not L<Bio::PrimarySeqI>)
objects can therefore be used for the nodes but any object that supports 
annotation objects and the object_id() method should work fine. 

At present it is fairly 'lightweight' in that it represents nodes and
edges but does not contain all the data about experiment ids etc. found
in the Protein Standards Initiative schema. Hopefully that will be
available soon.

A dataset may contain duplicate or redundant interactions. 
Duplicate interactions are interactions that occur twice in the dataset 
but with a different interaction ID, perhaps from a different 
experiment. The dup_edges method will retrieve these.

Redundant interaction are interactions that occur twice or more in a 
dataset with the same interaction id. These are more likely to be 
due to database errors. These methods are useful when merging 2 
datasets using the union() method. Interactions present in both 
datasets, with different IDs, will be duplicate edges. 

=head2 For Developers

In this module, nodes are represented by L<Bio::Seq::RichSeq> objects
containing all possible database identifiers but no sequence, as
parsed from the interaction files. However, a node represented by a
L<Bio::PrimarySeq> object should work fine too.

Edges are represented by L<Bio::Graph::Edge> objects. In order to
work with SimpleGraph these objects must be array references, with the
first 2 elements being references to the 2 nodes. More data can be
added in $e[2]. etc. Edges should  be L<Bio::Graph::Edge> objects, which 
are L<Bio::IdentifiableI> implementing objects.

At present edges only have an identifier and a weight() method, to 
hold confidence data, but subclasses of this could hold all the 
interaction data held in an XML document.

So, a graph has the following data:

1. A hash of nodes ('_nodes'), where keys are the text representation of a 
nodes memory address and values are the sequence object references.

2. A hash of neighbors ('_neighbors'), where keys are the text representation of a 
nodes memory address and a value is a reference to a list of 
neighboring node references.

3. A hash of edges ('_edges'), where a key is a text representation of the 2 nodes.
E.g., "address1,address2" as a string, and values are Bio::Graph::Edge 
objects.

4. Look up hash ('_id_map') for finding a node by any of its ids. 

5. Look up hash for edges ('_edge_id_map') for retrieving an edge 
object  from its identifier.

6. Hash ('_components').

7. An array of duplicate edges ('_dup_edges').

8. Hash ('_is_connected').

=head1  REQUIREMENTS

To use this code you will need the Clone.pm module availabe from CPAN.
You also need Class::AutoClass, available from CPAN as well.  To read in
XML data you will need XML::Twig available from CPAN.

=head1 SEE ALSO

L<Bio::Graph::SimpleGraph>
L<Bio::Graph::IO>
L<Bio::Graph::Edge>
L<Bio::Graph::IO::dip>
L<Bio::Graph::IO::psi_xml>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS

 Richard Adams - this module, Graph::IO modules.

 Email richard.adams@ed.ac.uk

=head2 AUTHOR2

 Nat Goodman - SimpleGraph.pm, and all underlying graph algorithms.

=cut

use strict;
package Bio::Graph::ProteinGraph;
use Bio::Graph::Edge;
use Clone qw(clone);
use base qw(Bio::Graph::SimpleGraph);

=head2  has_node

 name      : has_node
 purpose   : Is a protein in the graph?
 usage     : if ($g->has_node('NP_23456')) {....}
 returns   : 1 if true, 0 if false
 arguments : A sequence identifier.

=cut

sub has_node {

 my ($self, $arg) = @_;
 if (!$arg) {
		$self->throw ("I need a sequence identifier!");
   }
 my @nodes = $self->nodes_by_id($arg);
 if (defined($nodes[0])){return 1;}else{return 0};

}

=head2 nodes_by_id

 Name      : nodes_by_id
 Purpose   : get node memory address from an id
 Usage     : my @neighbors= $self->neighbors($self->nodes_by_id('O232322'))
 Returns   : a SimpleGraph node representation ( a text representation
             of a node needed for other graph methods e.g.,
             neighbors(), edges()
 Arguments : a protein identifier., e.g., its accession number.

=cut

sub nodes_by_id {

	my $self  = shift;
	my @nodes  = $self->_ids(@_);
	wantarray? @nodes: $nodes[0];

}


=head2 union

 Name        : union
 Purpose     : To merge two graphs together, flagging interactions as 
               duplicate.
 Usage       : $g1->union($g2), where g1 and g2 are 2 graph objects. 
 Returns     : void, $g1 is modified
 Arguments   : A Graph object of the same class as the calling object. 
 Description : This method merges 2 graphs. The calling graph is modified, 
               the parameter graph ($g2) in usage) is unchanged. To take 
               account of differing IDs identifying the same protein, all 
               ids are compared. The following rules are used to modify $g1.

               First of all both graphs are scanned for nodes that share 
               an id in common. 

         1. If 2 nodes(proteins) share an interaction in both graphs,
            the edge in graph 2 is copied to graph 1 and added as a
            duplicate edge to graph 1,

         2. If 2 nodes interact in $g2 but not $g1, but both nodes exist
            in $g1, the attributes of the interaction in $g2 are 
            used to make a new edge in $g1.

         3. If 2 nodes interact in g2 but not g1, and 1 of them is a new
            protein, that protein is put in $g1 and a new edge made to
            it. 

         4. At present, if there is an interaction in $g2 composed of a
            pair of interactors that are not present in $g1, they are 
            not copied to $g1. This is rather conservative but prevents
            the problem of having redundant nodes in $g1 due to the same
            protein being identified by different ids in the same graph.

         So, for example 

              Edge   N1  N2 Comment

    Graph 1:  E1     P1  P2
              E2     P3  P4
              E3     P1  P4

    Graph 2:  X1     P1  P2 - will be added as duplicate to Graph1
              X2     P1  X4 - X4 added to Graph 1 and new edge made
              X3     P2  P3 - new edge links existing proteins in G1
              X4     Z4  Z5 - not added to Graph1. Are these different
                              proteins or synonyms for proteins in G1?

=cut

sub union {

	my ($self, $other) = @_;
	my $class      = ref($self);
	if (!$other->isa($class)) {
		$self->throw("I need a ". $class . " object, not a [".
						 ref($other). "] object");
	}
	my @common_nodes;
	my %detected_common_nodes;
	my %seen_ids; # holds ids of nodes  already known to be common. 

	## for each node see if Ids are in common between the 2 graphs
	## just get1 common id per sequence.

	## Produces too many common nodes we only need 1 common id between nodes.
	for my $id (sort keys %{$self->{'_id_map'}}) {
		if (exists($other->{'_id_map'}{$id}) ) { 
			## check  if this node has a commonlink kown lready:
			my $node = $self->nodes_by_id($id);
			my $acc = $node->object_id;
			if (!exists($detected_common_nodes{$acc})) {
			   push @common_nodes, $id; ## we store the common id
			   $detected_common_nodes{$acc} = undef; ## this means we won't store >1 common identifier
			}
		}
	}

	## now cyle through common nodes..
	$self->debug( "there are ". scalar @common_nodes. " common nodes\n");
	my $i = 0;
	for my $common (@common_nodes) {
		if ($i++ % 10 ==0 ) {
			$self->debug(".");
		}
		## get neighbours of common node for self and other
		my @self_ns   = $self->neighbors($self->nodes_by_id($common));
		my @other_ns  = $other->neighbors($other->nodes_by_id($common));

		## now get all ids of all neighbours
		my %self_n_ids = $self->_get_ids(@self_ns);	# get all ids of neighbors

		## cycle through other neighbors
		for my $other_n(@other_ns){ 
			my %other_n_ids = $self->_get_ids($other_n); # get ids of single other neighbor

			## case (1) in description
			## do any ids in other graph exist in self ?
			# if yes,  @int_match is defined, interaction does not involve a new node
			my @int_match = grep{exists($self->{'_id_map'}{$_}) } keys %other_n_ids;
			if (@int_match){
				my $i = 0;
				my $edge;

				## we cycle through until we have an edge defined, this deals with 
				## multiple id matches
				while (!$edge && $i <= $#int_match){

					## get edge from other graph
					my $other_edge = $other->edge(
												 [$other->nodes_by_id($common),
												  $other->nodes_by_id($other_n->object_id)]
														  );

					## copy it
					my $edge = Bio::Graph::Edge->new(
													 -weight=> $other_edge->weight(),
													 -id    => $other_edge->object_id(),
													 -nodes =>[$self->nodes_by_id($common),
								   $self->nodes_by_id($int_match[$i])
																 ]);
					## add it to self graph.
					## add_edge() works out if the edge is a new,  
					## duplicate or a redundant edge.
					$self->add_edge($edge);

					$i++;
				}
			}		# end if
			## but if other neighbour is entirely new, clone it and 
			## make connection.
			else  {
				my $other_edge = $other->edge($other->nodes_by_id($other_n->object_id()),
														$other->nodes_by_id($common));
				my $new = clone($other_n);
				$self->add_edge(Bio::Graph::Edge->new(
							  -weight => $other_edge->weight(),
							  -id     => $other_edge->object_id(),
							  -nodes  =>[$new, $self->nodes_by_id($common)],
																 )
									);

				## add new ids to self graph look up table
				map {$self->{'_id_map'}{$_} = $new} keys %other_n_ids;
			} #end if
		} #next neighbor
	} #next node
}

=head2 edge_count

 Name     : edge_count
 Purpose  : returns number of unique interactions, excluding 
            redundancies/duplicates
 Arguments: void
 Returns  : An integer
 Usage    : my $count  = $graph->edge_count;

=cut

sub edge_count {

    my $self = shift;
    return scalar keys %{$self->_edges};

}

=head2 node_count

 Name     : node_count
 Purpose  : returns number of nodes.
 Arguments: void
 Returns  : An integer
 Usage    : my $count = $graph->node_count;

=cut

sub node_count {

    my $self = shift;
    return scalar keys %{$self->_nodes};

}

=head2 neighbor_count

 Name      : neighbor_count
 Purpose   : returns number of neighbors of a given node
 Usage     : my $count = $gr->neighbor_count($node)
 Arguments : a node object
 Returns   : an integer

=cut

sub neighbor_count{

 my ($self, $node) = @_;
 if (!$node->isa('Bio::NodeI')) {
  $self->throw ("I need a Bio::NodeI implementing object here , not a " . ref($node) . ".");
	}
 my @nbors = $self->neighbors($node);
 return scalar @nbors;
}

=head2 _get_ids_by_db

 Name     : _get_ids_by_db
 Purpose  : gets all ids for a node, assuming its Bio::Seq object
 Arguments: A Bio::SeqI object
 Returns  : A hash: Keys are db ids, values are accessions
 Usage    : my %ids = $gr->_get_ids_by_db($seqobj);

=cut

sub _get_ids_by_db {
	my %ids;
	my $dummy_self = shift;
	while (my $n = shift @_ ){  #ref to node, assume is a Bio::Seq
		if (!$n->isa('Bio::AnnotatableI') || ! $n->isa('Bio::IdentifiableI' )) {
			$n->throw("I need a Bio::AnnotatableI and Bio::IdentifiableI  implementing object, not a [" .ref($n) ."]");
		}

		##if BioSeq getdbxref ids as well.
		my $ac = $n->annotation();	
		for my $an($ac->get_Annotations('dblink')) {
			$ids{$an->database()} = $an->primary_id();
		}
	}
	return %ids;
}

sub _get_ids {

	my %ids;
	my $dummy_self = shift;
	while (my $n = shift @_ ){  #ref to node, assume is a Bio::Seq
		if (!$n->isa('Bio::AnnotatableI') || ! $n->isa('Bio::IdentifiableI' )) {
			$n->throw("I need a Bio::AnnotatableI and Bio::IdentifiableI  implementing object, not a [" .ref($n) ."]");
		}
		#get ids
		map {$ids{$_} = undef} ($n->object_id);

		##if BioSeq getdbxref ids as well.
		if ($n->can('annotation')) {
			my $ac = $n->annotation();	
			for my $an($ac->get_Annotations('dblink')) {
				$ids{$an->primary_id()} = undef;
			}
		}
	}
	return %ids;

}

=head2 add_edge

 Name        : add_edge
 Purpose     : adds an interaction to a graph.
 Usage       : $gr->add_edge($edge)
 Arguments   : a Bio::Graph::Edge object, or a reference to a 2 element list. 
 Returns     : void
 Description : This is the method to use to add an interaction to a graph. 
               It contains the logic used to determine if a graph is a 
               new edge, a duplicate (an existing interaction with a 
               different edge id) or a redundant edge (same interaction, 
               same edge id).

=cut

sub add_edge {

	my $self      = shift;
	my $edges     = $self->_edges;
	my $neighbors = $self->_neighbors;
	my $dup_edges = $self->_dup_edges;
	my $edge;
	while (@_) {
		if ( ref($_[0]) eq 'ARRAY' || !ref($_[0])) {
      	$self->SUPER::add_edges(@_);
			return;
		} 
		elsif ( $_[0]->isa('Bio::Graph::Edge') ) {	# it's already an edge
			$edge = shift;
		}
		else {
			$self->throw(" Invalid edge! - must be an array of nodes, or an edge object");
		}

		my ($m, $n) = $edge->nodes();
		next if $m eq $n;		# no self edges
		last unless defined $m && defined $n;
		($m,$n) = ($n,$m) if "$n" lt "$m";

		if (!exists($edges->{$m,$n})) {
			$self->add_node($m,$n);
			($m,$n)         = $self->nodes($m,$n);
			$edges->{$m,$n} = $edge;
			push(@{$neighbors->{$m}},$n);
			push(@{$neighbors->{$n}},$m);

			## create look up hash for edge ##
			$self->{'_edge_id_map'}{$edge->object_id()} = $edge;
		} else {
			## is it a redundant edge, ie with same edge id?
			my $curr_edge = $edges->{$m,$n};
			if($curr_edge->object_id() eq $edge->object_id()) {
				$self->redundant_edge($edge);
			}
			## else it is a duplicate i.e., same nodes but different edge id
			else {
				$self->add_dup_edge($edge); 
			}
		}
	}
	$self->_is_connected(undef);	# clear cached value

}

=head2 subgraph

 Name      : subgraph
 Purpose   : To construct a subgraph of  nodes from the main network.This 
             method overrides that of Bio::Graph::SimpleGraph in its dealings with 
             Edge objects. 
 Usage     : my $sg = $gr->subgraph(@nodes).
 Returns   : A subgraph of the same class as the original graph. Edge objects are 
             cloned from the original graph but node objects are shared, so beware if you 
             start deleting nodes from the parent graph whilst operating on subgraph nodes. 
 Arguments : A list of node objects.

=cut

sub subgraph {
 my $self=shift;

  ## make new graph of same type as parent
  my $class    = ref($self);
  my $subgraph = new $class;
  $subgraph->add_node(@_);
  # add all edges amongst the nodes
  my @nodes=$subgraph->nodes;
  my $i = 1;
  while(@nodes) {
    if ($i++ % 100 == 0) { print STDERR ".";}
    my $m=shift @nodes;
    my $edges = $self->_edges;
    for my $n (@nodes) { 
       if ($self->has_edge([$m,$n])) {
          my ($edge) = $self->edges([$m,$n]); ## returns list of edges
          my $id = $edge->object_id;
          $subgraph->add_edge(Bio::Graph::Edge->new(-nodes=>[$m,$n],
                                                    -id   => $id));
        }
    }
  }#next node
  return $subgraph;
}

=head2 add_dup_edge

 Name       : add_dup_edge
 Purpose    : to flag an interaction as a duplicate, take advantage of 
              edge ids. The idea is that interactions from 2 sources with 
              different interaction ids can be used to provide more 
              evidence for a interaction being true, while preventing 
              redundancy of the same interaction being present more than 
              once in the same dataset. 
 Returns    : 1 on successful addition, 0 on there being an existing 
              duplicate. 
 Usage      : $gr->add_dup_edge(edge->new (-nodes => [$n1, $n2],
                                           -score => $score
                                           -id    => $id);
 Arguments  : an EdgeI implementing object.
 Descripton : 


=cut

sub add_dup_edge {

	## get the 2 nodes
	my ($self, $newedge) = @_;
	## prelimaries
	my $newedge_id   = $newedge->object_id();

	## now we have node objects, an edge id.
	## is edge id new?
	my $dup_edges = $self->_dup_edges();
	if(!grep{$_->object_id eq $newedge_id } @$dup_edges) {
		push @$dup_edges, $newedge;
		}
	else {
		$self->redundant_edge($newedge);
	}
}

=head2 edge_by_id

 Name        : edge_by_id
 Purpose     : retrieve data about an edge from its id
 Arguments   : a text identifier
 Returns     : a Bio::Graph::Edge object or undef
 Usage       : my $edge = $gr->edge_by_id('1000E');

=cut

sub edge_by_id {

 my ($self, $id) = @_;
 if (!$id) {
	$self->warn ("Need a text identifier");
   	return;
	}
 if (ref($id)) {
    $self->throw(" I need a text identifier, not a [" . ref($id) . "].");
    }
  if (defined($self->{'_edge_id_map'}{$id})) {
     return $self->{'_edge_id_map'}{$id};
       }else {return;}

}


=head2 remove_dup_edges 

 Name        : remove_dup_edges
 Purpose     : removes duplicate edges from graph
 Arguments   : none         - removes all duplicate edges
               edge id list - removes specified edges
 Returns     : void
 Usage       :    $gr->remove_dup_edges()
               or $gr->remove_dup_edges($edgeid1, $edgeid2);

=cut

sub  remove_dup_edges{
  my ($self, @args) = @_;
  my $dups = $self->_dup_edges(); 
	if (!@args) {
  		$dups   = [];
		}
	else {
		while (my $node = shift @args) {
			my @new_dups;
			for my $dup (@$dups) {
				if (!grep{$node eq $_} $dup->nodes) {
					push @new_dups, $dup;
				}
			}
			$dups = \@new_dups;
		}
	}
	return 1;

}

=head2 redundant_edge

 Name        : redundant_edge
 Purpose     : adds/retrieves redundant edges to graph
 Usage       : $gr->redundant_edge($edge)
 Arguments   : none (getter) or a Biuo::Graph::Edge object (setter). 
 Description : redundant edges are edges in a graph that have the 
               same edge id, ie. are 2 identical interactions. 
               With edge arg adds it to list, else returns list as reference. 

=cut

sub redundant_edge {

	my ($self, $edge) =@_;
	if ($edge) {
		if (!$edge->isa('Bio::Graph::Edge')) {
			$self->throw ("I need a Bio::Graph::Edge object , not a [". ref($edge). "] object.");
		}
		if (!exists($self->{'_redundant_edges'})) {
			$self->{'_redundant_edges'} = [];
		}
		## add edge to list if not already listed
		if (!grep{$_->object_id eq $edge->object_id} @{$self->{'_redundant_edges'}}){
			push @{$self->{'_redundant_edges'}}, $edge;
		}
	}
	else {
		if (exists ($self->{'_redundant_edges'})){
			return @{$self->{'_redundant_edges'}};
		}
		else{
			
		}
	}
}

=head2 redundant_edges

 Name         : redundant_edges
 Purpose      : alias for redundant_edge

=cut

sub redundant_edges {
	my $self = shift;
	return $self->redundant_edge(shift);
}

=head2 remove_redundant_edges 

 Name        : remove_redundant_edges
 Purpose     : removes redundant_edges from graph, used by remove_node(),
               may be better as an internal method??
 Arguments   : none         - removes all redundant edges
               edge id list - removes specified edges
 Returns     : void
 Usage       :    $gr->remove_redundant_edges()
               or $gr->remove_redundant_edges($edgeid1, $edgeid2);

=cut

sub remove_redundant_edges {
my ($self, @args) = @_;
  my @dups = $self->redundant_edge(); 
	## if no args remove all 
	if (!@args) {
		$self->{'_redundant_edges'} = [];
		}
	else {
		while (my $node = shift @args) {
			my @new_dups;
			for my $dup (@dups) {
				if (!grep{$node eq $_} $dup->nodes) {
					push @new_dups, $dup;
				}
			}
			$self->{'_redundant_edges'} = \@new_dups;
		}
	}
	return 1;

}

=head2 clustering_coefficient

 Name      : clustering_coefficient
 Purpose   : determines the clustering coefficient of a node, a number 
             in range 0-1 indicating the extent to which the neighbors of
             a node are interconnnected.
 Arguments : A sequence object (preferred) or a text identifier
 Returns   : The clustering coefficient. 0 is a valid result.
             If the CC is not calculable ( if the node has <2 neighbors),
                returns -1.
 Usage     : my $node = $gr->nodes_by_id('P12345');
             my $cc   = $gr->clustering_coefficient($node);

=cut

sub clustering_coefficient {
	my  ($self, $val)  = @_;
	my $n = $self->_check_args($val);
	$self->throw("[$val] is an incorrect parameter, not presnt in the graph")
		unless defined($n);
	my @n = $self->neighbors($n);
	my $n_count = scalar @n;
	my $c = 0;

	## calculate cc if we can
	if ($n_count >= 2){
		for (my $i = 0; $i <= $#n; $i++ ) {
			for (my $j = $i+1; $j <= $#n; $j++) {
				if ($self->has_edge($n[$i], $n[$j])){
					$c++;
				}
			}
		}
		$c = 2 * $c / ($n_count *($n_count - 1));
		return $c; # can be 0 if unconnected. 
	}else{
		return -1; # if value is not calculable
	}
}

=head2 remove_nodes

 Name      : remove_nodes
 Purpose   : to delete a node from a graph, e.g., to simulate effect 
             of mutation
 Usage     : $gr->remove_nodes($seqobj);
 Arguments : a single $seqobj or list of seq objects (nodes)
 Returns   : 1 on success

=cut


sub remove_nodes {
	my $self = shift @_;
	if (!@_) {
		$self->warn("You have to specify a node");
		return;
		}
	my $edges     = $self->_edges;
	my $ns = $self->_neighbors;
	my $dups      = $self->_dup_edges;
	my $nodes     = $self->_nodes;
	while (my $val = shift @_ ) {

		## check argument
		my $node = $self->_check_args($val);
		$self->throw("[$val] is an incorrect parameter, not present in the graph")
				unless defined($node);

		##1. remove dup edges and redundant edges containing the node ##
		$self->remove_dup_edges($node);
		$self->remove_redundant_edges($node);

		##2. remove node from interactor's neighbours
		my @ns = $self->neighbors($node);
		for my $n (@ns) {
			my @otherns    = $self->neighbors($n); #get neighbors of neighbors 
			my @new_others = ();
			##look for node in neighbor's neighbors
			@new_others    = grep{$node ne $_} @otherns;
			$ns->{$n}   = \@new_others;
		}

		##3. Delete node from neighbour hash
		delete $ns->{$node};

		##4. Now remove edges involving node
		for my $k (keys %$edges) {
			##access via internal hash rather than by object. 
			if ($edges->{$k}->[0] eq $node ||
			   		$edges->{$k}->[1] eq $node){
                ## delete edge from look up hash
                my $edge_id = $edges->{$k}->object_id();
                delete $self->{'_edge_id_map'}{$edge_id};
				delete($edges->{$k});
			}
		}

		##5. Now remove node itself;
		delete $nodes->{$node}{'_node_id'};
		delete $nodes->{$node};

		##6. now remove aliases from look up hash so it can no longer be accessed.
		## is this wise? or shall we keep the sequence object available??
	}
	return 1;
}

=head2 unconnected_nodes

 Name      : unconnected_nodes
 Purpose   : return a list of nodes with no connections. 
 Arguments : none
 Returns   : an array or array reference of unconnected nodes
 Usage     : my @ucnodes = $gr->unconnected_nodes();

=cut

sub unconnected_nodes {
 my $self = shift;
 my $neighbours = $self->_neighbors;
 my $nodes      = $self->_nodes;
 my $uc_nodes   = [];
 for my $n (keys %$neighbours) {
	 if (@{$neighbours->{$n}} == 0){ 
		 push @$uc_nodes, $nodes->{$n};
	 }
 }
 wantarray?@$uc_nodes:$uc_nodes;
}

=head2 articulation_points

 Name      : articulation_points
 Purpose   : to find edges in a graph that if broken will fragment
               the graph into islands.
 Usage     : my $edgeref = $gr->articulation_points();
             for my $e (keys %$edgeref) {
				   print $e->[0]->accession_number. "-".
                     $e->[1]->accession_number ."\n";
             }
 Arguments : none
 Returns   : a list references to nodes that will fragment the graph 
             if deleted. 
 Notes     : This is a "slow but sure" method that works with graphs
               up to a few hundred nodes reasonably fast.

=cut

sub articulation_points {

 my $self      = shift;
 ## see if results are cahced already
 $self->{'_artic_points'} ||= '';
 return $self->{'_artic_points'} if $self->{'_artic_points'};

## else calculate...
 $self->debug( "doing subgraphs\n");
 my @subgraphs = $self->components();
 
 my %rts;

 for my $sg (@subgraphs) {
     my $all_nodes = $sg->_nodes;
     $self->debug( "in subgraph - size". scalar (keys %$all_nodes) . "\n");
     ##ignore isolated vertices
     next if scalar keys %$all_nodes <= 2;
     my $neighbors = $sg->_neighbors;

     ## find most connected - will be artic point if has >2 neighbors.
     ## use this to initiate DFS
     my ($c, $id);
     my $max = 0;
     for my $n (keys %$neighbors) {
	 my $c = scalar @{$neighbors->{$n}};#
	 ($max, $id) = ($c, $n) if  $c > $max;#
     }

     my $t      = $sg->node_traversal($all_nodes->{$id},'d');
     my @nodes  = $t->get_all();
     $id = 0;
     #assign node ids
     for my $n(@nodes) {
	 $n->{'_node_id'} = $id;	
	 $id++;
     }

     ## cycle through each node 
     for (my $i       = $#nodes; $i >= 0; $i--) {

	 ## initiate minimumn to node_id
	 my $curr_min = $all_nodes->{$nodes[$i]}{'_node_id'};
	 #print STDERR "currmin - $curr_min, i = $i\n";
	 ## cycle through neighbors, reset minumum if required
	 my $nbors    = $neighbors->{$nodes[$i]};
	 for my $nbor (@$nbors) {	
	     my $nbor_id = $all_nodes->{$nbor}{'_node_id'};

	     ## if is back edge ##
	     if ($nbor_id < $i) {
		 $curr_min  = $nbor_id if $nbor_id < $curr_min ;
	     }

	     ## else is tree edge
	     elsif($nbor_id > $i) {
		 my $wlow   = $all_nodes->{$nbor}{'_wlow'};
		 $curr_min  = $wlow if $wlow < $curr_min;
	     }
	 }#next neighbor

	 ## now we know the minimum, save. 
	 $all_nodes->{$nodes[$i]}{'_wlow'} = $curr_min;

	 ## now get tree nodes and test condition
	 my @treenodes = grep{$all_nodes->{$_}{'_node_id'} > $i}@$nbors;
	 for my $tn (@treenodes) {
	     if(($all_nodes->{$tn}{'_wlow'} >= $i && $i != 0) ||
		($i == 0  && scalar @{$neighbors->{$nodes[0]}} > 1) ){
		 $rts{$nodes[$i]} = $nodes[$i] unless exists $rts{$nodes[$i]};
	     }
	 }

     }#next node
 }#next sg
## cache results and return
$self->{'_artic_points'} =   [values %rts]; ## 
return $self->{'_artic_points'}; 
}

=head2 is_articulation_point

 Name      : is_articulation_point
 Purpose   : to determine if a given node is an articulation point or not. 
 Usage     : if ($gr->is_articulation_point($node)) {.... 
 Arguments : a text identifier for the protein or the node itself
 Returns   : 1 if node is an articulation point, 0 if it is not 

=cut

sub is_articulation_point {
	my ($self, $val) = @_;
	my $node = $self->_check_args($val);
 
	## this uses a cached value so it does not have to recalculate each time..
	my $artic_pt_ref = $self->articulation_points();
	my $acc = $node->accession_number;
	if (grep{$_->accession_number eq $acc} @$artic_pt_ref ){
		return 1;
   }
	else {
		return 0;
   }
}

sub _ids {
	my $self = shift;
	my @refs;
	while (my $id = shift@_) {
		push @refs, $self->{'_id_map'}{$id};
	}
	return @refs;
}

sub _check_args {
## used to check a parameter is a valid node or a text identifier
	my ($self, $val) = @_;
	my $n;
	if (!$val ) {
		$self->throw( "I need a node that's a Bio::AnnotatableI and Bio::IdentifiableI");
		}

	## if param is text try to get sequence object..
	if (!ref($val)){
		 $n = $self->nodes_by_id($val);
		if(!defined($n)) {
			$self->throw ("Cannnot find node given by the id [$val]");
			}
	}
	# if reference should be a NodeI implementing object.
    elsif (!$val->isa('Bio::AnnotatableI') || !$val->isa('Bio::IdentifiableI')) {
		$self->throw( "I need a node that's a Bio::AnnotatableI and Bio::IdentifiableI ,not a [". ref($val) . "].");
		}

	## is a seq obj
	else {$n = $val};
	return $n; #n is either a node or undef
}

1;
