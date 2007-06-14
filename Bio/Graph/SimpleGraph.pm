# $Id$
#
# BioPerl module for Bio::Graph::SimpleGraph
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::SimpleGraph - create and manipulate undirected graphs

=head1 SYNOPSIS

  use Bio::Graph::SimpleGraph;

  my $graph=new SimpleGraph;
  # read pairs of nodes from STDIN
  while (<>) {
    my($node1,$node2)=split;
    $graph->add_edge($node1,$node2);
  }
  my @nodes=graph->nodes;	    # get list of nodes
  my @edges=graph->edges;	    # get list of edges
  foreach my $node (@nodes) {
    my @neighbors=$node->neighbors; # get list of neighboring nodes
  }

=head1 DESCRIPTION

This is a simple, hopefully fast undirected graph package. The only
reason this exists is that the standard CPAN Graph pacakge,
Graph::Base, is seriously broken.  The package implements a small and
eclectic assortment of standard graph algorithms that we happened to
need for our applications.

This module is a subclass of Class::AutoClass (available at CPAN).
AutoClass auotgenerates simple accessor and mutator methods (aka get
and set methods).  It also automates class initialization.

Nodes can be any Perl values, including object references. Edges are 
pairs of nodes. 

(Caveat: be careful with values that contain embedded instances of $;
(the character Perl uses to separate components of multi-dimensional
subscripts), because we use this in the text representation of edges.

The main data structures are:

  An edge (x,y) is represented canonically as a two element list in
  which the lexically smaller value is first.  Eg, the node ('b','a')
  is represented as ['a','b'].  

  The graph contains 

  1) A hash mapping the text representation of a node to the node
     itself.  This is mostly relevant when the node is a reference.

  2) A hash mapping the text representation of a node to a list of 
     the node's neighbors.

  3) A hash mapping the text representation of an edge to the edge itself.


=head1 KNOWN BUGS AND CAVEATS

This is still a work in progress.

=head1 AUTHOR - Nat Goodman

Email natg@shore.net

=head1 COPYRIGHT

Copyright (c) 2003 Institute for Systems Biology (ISB). All Rights Reserved.
This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=head1 APPENDIX

=head2 Conventions for nodes and edges

A node can be any Perl values, including an object, ARRAY, or HASH
reference.  When nodes are references, the software often works with
the text representaion of the reference, ie, what you get if you print
the reference.  This can be confusing.  Sorry.  For example if a node
is the HASH

  {name=>'caspase-9',symbol=>'CASP9'}

The text representation would be something like

  HASH(0x804c830)

When nodes are scalar values, eg, a string, the value and the text
representation are the same.  This is a common case in test programs
and examples, but less common in real applications.

An edge is represented internally as an ARRAY ref of two nodes, in
which the lexically smaller value is first.  Actually, the first node
is the one whose text representation is lexically smaller.

When passing edges as arguments to SimpleGraph methods, the edge can
be represented in several ways.

  1) An ARRAY ref of the nodes, eg, ['a','b'].

  2) A list of the two nodes, eg, ('a','b')

  3) Form (1) or (2) using the text represention of the node 
     instead of the node itself

You needn't worry about which node is lexically smaller.  SimpleGraph
performs this calculation internally.

When SimpleGraph returns edges as results, they are always in form
(1), ie, as ARRAY refs of nodes in correct lexical order.

=head2 General conventions for methods

When methods return lists, we generally check the context (via
wantarray) and return an ARRSY or ARRAY ref as appropriate.  We're not
100% consistent in this (sorry), so check the code if you have doubts.

We often define singular and plural forms of methods, eg, node and
nodes.  These differ in how they behave in a scalar context.  The
singular form assumes you want one answer and returns that, while the
plural form assumes you want a list of answers are returns it as an
ARRAY ref.  We're not 100% consistent in this (sorry), so check the
code if you have doubts.

The rest of the documentation describes the methods.

=head2 Constructors

 Title   : new (inherited from Class::AutoClass)
 Usage   : my $graph=new SimpleGraph;
 Function: Create new SimpleGraph object
 Returns : Newly created object
 Args    : (optional)
           nodes=>ARRAY of nodes, eg, ['a','b','c']
           edges=>ARRAY of edges, see add_edges for details

=head2 Basic node and edge operations

 Title   : add_nodes, add_node
 Usage   : $graph->add_nodes('a','b');
           $graph->add_node('a')) {
 Function: Add nodes to graph. Nodes that are already in graph
           are ignored.
 Args    : ARRAY of nodes.
 Returns : Nothing useful
 Note    : Singular and plural forms are synonymous

 Title   : add_edges, add_edge
 Usage   : $graph->add_edges('a','b',['b','c']);
           $graph->add_edge('c','d')) {
 Function: Add edges to graph. 
           Edges that are already in graph are not added again, but
           are placed in a separate 'duplicate edges' list.
           Automatically adds any nodes that are not yet in the graph.
 Args    : ARRAY of edges in any of the forms described in the
           previous section.  The forms can be mixed as shown in
           the Usage here.
 Returns : Nothing useful
 Note    : Singular and plural forms are synonymous

 Title   : nodes
 Usage   : my @nodes=$graph->nodes;
           if (@{$graph->nodes('a','b')}==2) {
             print "a, b are both nodes\n";
           }
 Function: Return all nodes or the given ones.  
           With no args returns all nodes.  
           With args, returns the nodes corresponding to each arg, or
           undef if the arg is not a node.  Useful for testing whether
           a given value is a node in the graph.
 Args    : (optional)
           ARRAY of nodes or text representations of nodes
 Returns : ARRAY or ARRAY ref of nodes (for args that correspond to
           nodes), or undef (for args that are not nodes)

 Title   : edges
 Usage   : my @edges=$graph->edges;
           if (@{$graph->edges('a','b',['b','c'])}==2) {
             print "[a,b] and [b,c] are both edges\n";
           }
 Function: Return all edges or the given ones.  
           With no args returns all edges.  
           With args, returns the edges corresponding to each arg, or
           undef if the arg is not a edge.  Useful for testing whether
           a given value is a edge in the graph.
 Args    : (optional)
           One or more edges in any of the forms described in the
           previous section.  The forms can be mixed as shown in
           the Usage here.
 Returns : ARRAY or ARRAY ref of edges for args that correspond to
           edges), or undef (for args that are not edges)

 Title   : node
 Usage   : if ($graph->node('a')) {
             print "a is a node\n";
           }
 Function: Test whether a value is a node in the graph, or map the
           text representation of a node to the node itself.  The
           method can also be fed a list of values (like the 'nodes'
           method) and it will test all of them.
 Args    : Usually, a single node.
           The function also accepts a list of nodes.
 Returns : In scalar context (the usual case): the node corresponding
           to the arg (if there's just one), or the node corresponding
           to the first arg (if a list of args were provided, which is
           kind of dumb in this case), or undef if the arg is not a
           node.

           In array context, it behaves just like 'nodes', returning
           an ARRAY of nodes (for args that correspond to nodes), or
           undef (for args that are not nodes)

 Title   : edge
 Usage   : if ($graph->edge('a','b')) {
             print "a,b is a edge\n";
           }
           if ($graph->edge(['a','b'])) {
             print "[a,b] is a edge\n";
           }
 Function: Test whether a value is a edge in the graph, or map the
           text representation of a edge to the edge itself.  The
           method can also be fed a list of edges (like the 'edges'
           method) and it will test all of them.
 Args    : Usually, a single edge.  Same format as 'edges'
           The function also accepts a list of edges, exactly like 
           'edges'
 Returns : In scalar context (the usual case): the edge corresponding
           to the arg (if there's just one), or the or the edge
           corresponding to the first arg (if a list of args were
           provided, which is kind of dumb in this case), or undef if
           the arg is not a edge.

           In array context, it behaves just like 'edge's, returning
           an ARRAY of edges (for args that correspond to edges), or
           undef (for args that are not edges)

 Title   : has_nodes, has_node
 Usage   : if ($graph->has_nodes('a','b')) {
             print "a, b are both nodes\n";
           }
           if ($graph->has_node('a')) {
             print "a is a node\n";
           }
 Function: Return true is all args are nodes.
 Args    : ARRAY of nodes or text representations of nodes
 Returns : Boolean
 Note    : Singular and plural forms are synonymous

 Title   : has_edges
 Usage   : if ($graph->has_edges('a','b',['b','c'])) {
             print "[a,b] and [b,c] are both edges\n";
           }
           if ($graph->has_edge('a','b')) {
             print "[a,b] is an edge\n";
           }
 Function: Return true is all args are edges.
 Args    : ARRAY of edges in the forms described in the section above
 Returns : Boolean
 Note    : Singular and plural forms are synonymous

 Title   : neighbors, neighbor
 Usage   : my @nodes=$graph->neighbors($node)
           my @nodes=$graph->neighbors($node,'node')
           my @edges=$graph->neighbors($edge,'edge');
 Function: Return the node or edge neighbors of a given node or edge.
 Args    : (mandatory)
           $source: node or edge whose neighbors are sought
           (optional)
           $what: the word 'node' or 'edge' (actually, anything starting
                  with 'n' or 'e' will do)
                  default: 'node'
 Returns : ARRAY or ARRAY ref of nodes or edges
 Note    : Singular and plural forms are synonymous. This may not be
           right.

 Title   : dup_edges
 Usage   : my @dups=$graph->dup_edges;
 Function: Return duplicate edges
 Args    : None
 Returns : ARRAY or ARRAY ref of edges that have been added more than
           once.

=head2 Graph properties

 Title   : is_connected
 Usage   : if ($graph->is_connected) {
             print "graph has only one connected component\n";
           }
 Function: Return true if the graph is connected
 Args    : None
 Returns : Boolean

 Title   : is_empty
 Usage   : if ($graph->is_empty) {
             print "graph has no nodes or edges\n";
           }
 Function: Return true if the graph is empty, ie, has no nodes or edges
 Args    : None
 Returns : Boolean

 Title   : is_tree
 Usage   : if ($graph->is_tree) {
             print "graph is a tree\n";
           }
 Function: Return true if the graph is a tree, ie, it's connected and
           has no cycles
 Args    : None
 Returns : Boolean

 Title   : is_forest
 Usage   : if ($graph->is_forest) {
             print "graph is a forest\n";
           }
 Function: Return true if the graph is a forest, ie, it has no cycles
           but may not be connected
 Args    : None
 Returns : Boolean

 Title   : is_cyclic
 Usage   : if ($graph->is_cyclic) {
             print "graph contains at least one cycle\n";
           }
 Function: Return true if the graph is a cyclic.
 Args    : None
 Returns : Boolean

 Title   : density
 Usage   : my $density=$graph->density
 Function: Compute graph 'density' which is the number of edges
           divided by the maximum possible number of edges
 Args    : None
 Returns : Number

=head2 Graph operations

 Title   : subgraph
 Usage   : my $subgraph=$graph->subgraph('a','b','c');
 Function: Compute node subgraph. Constructs a new graph whose nodes
           are the arguments, and whose edges are the edges of the
           original graph that only involve the given nodes.
 Args    : ARRAY of nodes or text representations of nodes
 Returns : New graph

 Title   : neighbor_subgraph
 Usage   : my $subgraph=$graph->subgraph('a');
 Function: Construct node subgraph graph whose nodes are the given
           node and its neighbors.  are the arguments, and whose edges
           are the edges of the original graph that only involve the
           given nodes.
 Args    : Node or text representations of node
 Returns : New graph

 Title   : union
 Usage   : my $union=$graph->union($other_graph);
 Function: Construct new graph whose nodes are the union of the nodes
           of the current graph and $other_graph, and whose edges are
           the union of the edges of the current graph and
           $other_graph.
 Args    : $other_graph: a graph
 Returns : New graph

 Title   : intersection
 Usage   : my $intersection=$graph->intersection($other_graph);
 Function: Construct new graph whose nodes are the intersection of the
           nodes of the current graph and $other_graph, and whose
           edges are the intersection of the edges of the current
           graph and $other_graph.
 Args    : $other_graph: a graph
 Returns : New graph

=head2 Graph algorithms

 Title   : traversal
 Usage   : my $traversal=$graph->traversal('a','depth first','node');
           my @nodes;
           while (my $node=$traversal->get_next) {
             push(@nodes,$node);
           }
           my $traversal=$graph->traversal('a','depth first','node');
           my @nodes=$traversal->get_all;
 Function: Do node or edge traversal in depth or breadth first order.
 Args    : (optional)
           $start: starting node or edge for traversal
                   default: software picks arbitrary start
           $order: 'depth first' or 'breadth first' (actually,
                   anything starting with 'd' or 'b' will do)
                   default: 'depth first'
           $what: 'node' or 'edge' (actually, anything starting
                  with 'n' or 'e' will do)
                  default: 'node'
 Returns : SimpleGraph::Traversal object
           This is an iterator with the following methods:

           get_next: get next item in traversal or undef if 
                     traversal is exhausted
           get_this: get current item in traversal
           get_all : get all remaining items in traversal as
                     ARRAY (in array context) or ARRAY ref
           has_next: return true if there are more items in
                     traversal, else undef
           reset   : restart traversal

 Note    : It's also possible, and perhaps easier, to perform a
           traversal by creating a SimpleGraph::Traversal object
           directly.  The constructor is

           new SimpleGraph::Traversal(-graph=>$graph,-start=>$start,
                                      -order=>$order,-what=>$what)

 Title   : node_traversal
 Usage   : my $traversal=$graph->node_traversal('a','depth first');
           my @nodes;
           while (my $node=$traversal->get_next) {
             push(@nodes,$node);
           }
           my $traversal=$graph->node_traversal('a','depth first');;
           my @nodes;


           my @nodes;

           my @nodes=$traversal->get_all;
 Function: Do node traversal in depth or breadth first order.
           Wrapper for 'traversal' method. See above.
 Args    : (optional)
           $start: starting node for traversal
                   default: software picks arbitrary start
           $order: 'depth first' or 'breadth first' (actually,
                   anything starting with 'd' or 'b' will do)
                   default: 'depth first'
 Returns : SimpleGraph::Traversal object

 Title   : edge_traversal
 Usage   : my $traversal=$graph->edge_traversal('a','depth first');
           my @edges;
           while (my $edge=$traversal->get_next) {
             push(@edges,$edge);
           }
           my $traversal=$graph->edge_traversal('a','depth first');
           my @edges=$traversal->get_all;
 Function: Do edge traversal in depth or breadth first order.
           Wrapper for 'traversal' method. See above.
 Args    : (optional)
           $start: starting edge for traversal
                   default: software picks arbitrary start
           $order: 'depth first' or 'breadth first' (actually,
                   anything starting with 'd' or 'b' will do)
                   default: 'depth first'
 Returns : SimpleGraph::Traversal object

 Title   : components
 Usage   : my @components=$graph->components;
           for my $component (@components) {
             my @nodes=$component->nodes;
             my @edges=$component->edges;
           }
 Function: Compute the connected components of the graph.  A connected
           component is a maximal connected subgraph.  'Connected'
           means you can get from any node of the component to any
           other by following a path.  'Maximal' means that every node
           you can reach from the component is in the component.
 Args    : None
 Returns : ARRAY or ARRAY ref of SimpleGraphs
 Note    : The software caches the components once computed, so it's efficient
           to call this repeatedly.

 Title   : shortest_paths
 Usage   : my @paths=$graph->shortest_paths;
           for my $path (@paths) {
             my @nodes_on_path=@$path;
             my $start=$nodes_on_path[0];
             my $end=$nodes_on_path[$#nodes_on_path];
           }
 Function: Compute shortest path between each pair of nodes.
 Args    : None
 Returns : ARRAY or ARRAY ref of paths, where each path is an ARRAY
           ref of nodes.  The result contains one path for each pair
           of nodes for which a path exists.

 Title   : connected_nodesets
 Usage   : my @nodesets=$graph->connected_nodesets;
           for my $nodeset (@nodesets) {
             my @nodes=@$nodeset;
           }
 Function: Compute all sets of nodes that form connected subgraphs. 
           A connected nodeset is a set of nodes such that it's
           possible to get from any node to any other by following a
           path that only includes nodes in the set. 
 Args    : None
 Returns : ARRAY or ARRAY ref of nodeset, where each nodeset is an ARRAY
           ref of nodes.  
 Note    : Use with caution.  The number of nodesets is very
           large for graphs that are highly connected.

 Title   : connected_subgraphs
 Usage   : my @subgraphs=$graph->connected_subgraphs;
 Function: Compute all connected subgraphs of the current graph.
 Args    : None
 Returns : ARRAY or ARRAY ref of subgraphs
 Note    : Use with caution.  The number of connected subgraphs is
           very large for graphs that are highly connected.

=cut

package Bio::Graph::SimpleGraph;
use vars qw(@AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS %DEFAULTS);
use Bio::Graph::SimpleGraph::Traversal;
use strict;
use base qw(Class::AutoClass);
@AUTO_ATTRIBUTES=qw(_nodes _edges _neighbors _dup_edges _is_connected _components);
%SYNONYMS=();
@OTHER_ATTRIBUTES=qw();
%DEFAULTS=(_nodes=>{},_edges=>{},_neighbors=>{},_dup_edges=>[]);
Class::AutoClass::declare(__PACKAGE__);

# Implementation:
#  An edge (x,y) is represented canonically as a two element list in which the 
#    lexically smaller value is first.  Eg, the node ('b','a') is represented 
#    as ['a','b'].
#  Graph contains
#    A hash mapping the text representation of a node to the node itself
#    A hash mapping the text representation of a node to the nodes neighbors
#    A hash mapping the text representation of an edge to the edge itself.

sub _init_self {
  my($self,$class,$args)=@_;
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  my $nodes = $args->nodes;
  defined($nodes) and $self->add_nodes(@$nodes);
  my $edges = $args->edges;
  defined($edges) and $self->add_edges(@$edges);
}
sub nodes {
  my $self = shift;
  my @ret  = @_? @{$self->_nodes}{@_}: values %{$self->_nodes};
  wantarray? @ret: \@ret;
}

sub node {
  my $self   = shift; 
  my @result = $self->nodes(@_); 
  wantarray? @result: $result[0];
}

sub edges {
  my $self = shift;
  my @ret;
  unless (@_) {
    @ret = values %{$self->_edges};
  } else {
    my $edges=$self->_edges;
    while (@_) {
      my($m,$n);
      if ('ARRAY' eq ref $_[0] || $_[0]->isa('Bio::Graph::Edge')) {	# it's already an edge
			my $edge = shift;
			($m,$n)  = @$edge[0..1];
      	} else {
		($m,$n)=(shift,shift);
      }
      last unless defined $m && defined $n;
      ($m,$n) = ($n,$m) if "$n" lt "$m";
      push(@ret, $edges->{$m,$n});
    }
  }
  wantarray? @ret: \@ret;
}

sub edge {
  my $self   = shift; 
  my @result = $self->edges(@_);
  wantarray? @result: $result[0];
}

sub has_nodes {
  my $self = shift;
  my @ret  = @{$self->_nodes}{@_};
  return (grep {!defined $_} @ret)? undef: 1;
}
sub has_node {my $self=shift; $self->has_nodes(@_); }

sub has_edges {
  my $self  = shift;
  my $edges = $self->_edges;
  while (@_) {
    my($m,$n);
    if ('ARRAY' eq ref $_[0] || $_[0]->isa('Bio::Graph::Edge')) {	# it's already an edge
      my $edge = shift;
      ($m,$n)  = @$edge[0..1];#first 2 elements are nodes
    } else {
      ($m,$n)  = (shift,shift);
    }
    last unless defined $m && defined $n;
    ($m,$n)=($n,$m) if "$n" lt "$m";
    return unless ($edges->{$m,$n});
  }
  return 1;
}
sub has_edge {my $self=shift; $self->has_edges(@_); }

sub neighbors {
  my($self,$source,$what)=@_;
  $what or $what='node';
  my $result;
  if ($what=~/^n/i) {
    $result=$self->_neighbors->{$source};
  } elsif ($what=~/^e/i) {
    my $edge=$self->edge($source);
    return unless $edge;
    my($m,$n)=@$edge;
    my %edges;
    for my $node (@{$self->_neighbors->{$m}}) {
      my $edge=$self->edge([$m,$node]);
      $edges{$edge}=$edge;
    }
    for my $node (@{$self->_neighbors->{$n}}) {
      my $edge=$self->edge([$n,$node]);
      $edges{$edge}=$edge;
    }
    delete $edges{$edge};	# remove source edge from result
    @$result=values %edges;
  } else {
    $self->throw("Unrecognized \$what parameter $what: should be 'node' or 'edge'");
  }
 if ($result){
 	 return wantarray? @$result: $result;
	}
	else { return ();}
}
sub neighbor {my $self=shift; $self->neighbors(@_); }

sub add_nodes {
  my $self      = shift @_;
  my $nodes     = $self->_nodes;
  my $neighbors = $self->_neighbors;
  for my $n (@_) {
    next if defined $nodes->{$n};
    $nodes->{$n}     = $n;
    $neighbors->{$n} = [];
  }
  $self->_components(undef);	# clear cached value
}
sub add_node {my $self=shift; $self->add_nodes(@_); }

sub add_edges {
  my $self  = shift @_;
  my $edges = $self->_edges;
  my $neighbors = $self->_neighbors;
  my $dup_edges = $self->_dup_edges;
  while (@_) {
    my($m,$n);
    if ('ARRAY' eq ref $_[0] ) {	# it's already an edge
      my $edge = shift;
      ($m,$n)  = @$edge;
    } else {
      ($m, $n )=(shift,shift, shift);
    }
    next if $m eq $n;		# no self edges
    last unless defined $m && defined $n;
    ($m,$n)=($n,$m) if "$n" lt "$m";
    unless ($edges->{$m,$n}) {
      $self->add_node($m,$n);
      ($m,$n) = $self->nodes($m,$n);
      $edges->{$m,$n} = [$m,$n];
      push(@{$neighbors->{$m}},$n);
      push(@{$neighbors->{$n}},$m);
    } else {
      push(@$dup_edges,[$m,$n]);
    }
  }
  $self->_is_connected(undef);	# clear cached value
}
sub add_edge {my $self=shift; $self->add_edges(@_); }

sub dup_edges {
  my $self=shift;
  my @ret=@{$self->_dup_edges};
  wantarray? @ret: \@ret;
}

sub is_connected {@{$_[0]->components}<=1;} # connected graph has 0 or 1 components
sub is_empty {@{$_[0]->nodes}==0;}

sub is_tree {			# tree is (connected and #edges=#nodes-1) or empty
  ($_[0]->is_connected && (@{$_[0]->edges}==(@{$_[0]->nodes}-1))) ||
    $_[0]->is_empty;
}
sub is_forest {
  my($self)=@_;
  my @components=$self->components;
  for my  $component (@components) {
    return unless $component->is_tree;
  }
  return 1;
}
sub is_cyclic {!$_[0]->is_forest}

sub density {
  my($self) = @_;
  my @nodes = $self->nodes;
  my @edges = $self->edges;
  return 0 if @nodes < 2;
  my $max = @nodes*(@nodes-1)/2;
  @edges/$max;
}

sub subgraph {
  my $self=shift;

  ## make new graph of same type as parent
  my $class    = ref($self);
  my $subgraph = new $class;

  $subgraph->add_node(@_);
  # add all edges amongst the nodes
  my @nodes=$subgraph->nodes;
  while(@nodes) {
    my $m=shift @nodes;
    for my $n (@nodes) {
      $subgraph->add_edge([$m,$n]) if $self->has_edge([$m,$n]);
    }
  }
#  my @edges=grep{my($n1,$n2)=@$_; $subgraph->has_node($n1,$n2)} $self->edges;
#  $subgraph->add_edge(@edges);
  $subgraph;
}
sub neighbor_subgraph {
  my($self,$n)=@_;
  $self->subgraph($self->node($n),$self->neighbors($n));
}
sub union {
  my($self,$other)=@_;
  my $result=Bio::Graph::SimpleGraph->new();
  $result->add_node($self->nodes,$other->nodes);
  $result->add_edge($self->edges,$other->edges);
  $result;
}
sub intersection {
  my($self,$other)=@_;
  my $result=Bio::Graph::SimpleGraph->new();
  for my $node ($self->nodes) {
    next unless $other->has_node($node);
    $result->add_node($node);
  }
  for my $edge ($self->edges) {
    next unless $other->has_edge($edge);
    $result->add_edge($edge);
  }
  $result;
}

sub traversal {
  my($self,$start,$order,$what)=@_;
  Bio::Graph::SimpleGraph::Traversal->new(-graph=>$self,-start=>$start,-order=>$order,-what=>$what);
}
sub node_traversal {
  my($self,$start,$order)=@_;
  Bio::Graph::SimpleGraph::Traversal->new(-graph=>$self,-start=>$start,-order=>$order,-what=>'node');
}
sub edge_traversal {
  my($self,$start,$order)=@_;
  Bio::Graph::SimpleGraph::Traversal->new(-graph=>$self,-start=>$start,-order=>$order,-what=>'edge');
}

sub shortest_paths {
  my($self)=@_;
  # initialization
  my @nodes=$self->nodes;
  my $dist={};
  my $path={};
  for (my $i=0; $i<@nodes; $i++) { # start from i
    my $node0=$nodes[$i];
    for (my $j=0; $j<@nodes; $j++) { # end at j
      my $node1=$nodes[$j];
      if ($i==$j) {
	$dist->{$i,$j}=0;
	next;
      }
      next unless $self->has_edge([$node0,$node1]);
      $dist->{$i,$j}=1;
      $path->{$i,$j}=[$i,$j];
    }
  }
  # compute paths
  for (my $k=0; $k<@nodes; $k++) {        # k is intermediate point
    for (my $i=0; $i<@nodes-1; $i++) {    # start from i
      next unless defined $dist->{$i,$k};
      for (my $j=$i+1; $j<@nodes; $j++) { # NG 04-02-10 added optimization
#      for (my $j=0; $j<@nodes; $j++) { # end at j
	next unless defined $dist->{$k,$j};
	# path i..k..j exists -- is it shorter than what we already have?
	if (!defined $dist->{$i,$j} || $dist->{$i,$k}+$dist->{$k,$j} < $dist->{$i,$j}) {
	  $dist->{$i,$j}=$dist->{$i,$k}+$dist->{$k,$j};
	  $path->{$i,$j}=_sp_join_paths($path,$i,$k,$j);
#	  # NG 04-02-10 next two lines needed for optimization above
	  $dist->{$j,$i}=$dist->{$i,$j};
	  $path->{$j,$i}=[reverse @{$path->{$i,$j}}];
	}
      }
    }
  }
  # convert node indices (i,j,..) into nodes
  my $paths={};
  for (my $i=0; $i<@nodes-1; $i++) {    # start from i
    my $nodei=$nodes[$i];
    for (my $j=$i+1; $j<@nodes; $j++) { # end at j
      my $p=$path->{$i,$j};
      my $nodej=$nodes[$j];
      my $path_nodes=[map {$nodes[$_]} @$p];
      if ("$nodei" lt "$nodej") {
	$paths->{$nodei,$nodej}=$path_nodes;
      } else {
	$paths->{$nodej,$nodei}=[reverse @$path_nodes];
      }
    }
  }
  my @paths=grep {@$_} values %$paths;
  wantarray? @paths : \@paths;
}
sub _sp_join_paths {
  my($path,$i,$k,$j)=@_;
  my $path0=$path->{$i,$k} || [];
  my $path1=$path->{$k,$j} || [];
  my $last0=@$path0-1;
  my $last1=@$path1-1;
  my $result=[];
  @$result=(@$path0[0..$last0-1],$k,@$path1[1..$last1]);
  $path->{$i,$j}=$result;
}

sub connected_nodesets {
  my($self)=@_;
  my @nodes=$self->nodes;
  my %node2num;
  # algorithm works with node numbers, not nodes
  for (my $i=0;$i<@nodes;$i++) {
    $node2num{$nodes[$i]}=$i;
  }
  my %nodesets;			     # nodesets that have been seen (not necessarily processed)
  my @future=map {[$_]} 0..@nodes-1; # nodesets to be processed -- init to nodes
  @nodesets{0..@nodes-1}=@future;

  while (@future) {
    my $present=shift @future;
    my %present;
    @present{@$present}=(1)x@$present; # for quick check of redundant neighbors
    for my $node (@$present) {
      for my $neighbor (map {$node2num{$_}} $self->neighbors($nodes[$node])) {
	next if $present{$neighbor};
	my $future=[sort by_num (@$present,$neighbor)];
	my $key=join($;,@$future);
	next if $nodesets{$key}; # skip if already seen
	$nodesets{$key}=$future;
	push(@future,$future);
      }
    }
  }
  my $nodesets=[];
  for my $nodeset (values %nodesets) {
    push(@$nodesets,[map {$nodes[$_]} @$nodeset]);
  }
  wantarray? @$nodesets: $nodesets;
}

sub connected_subgraphs {
  my($self)=@_;
  my $subgraphs;
  @$subgraphs=map {$self->subgraph(@$_)} $self->connected_nodesets;
  wantarray? @$subgraphs: $subgraphs;
}

sub components {
    my($self)  = @_;
    return $self->_components if defined $self->_components;
    my $components = [];
    my @nodes      = $self->nodes;
    my %future;
    my $i = 1;
    @future{@nodes}=(0)x@nodes;
    while(my($node, $used)=each %future) {
	if ($i++ %10 ==0 ) {
	    $self->debug("|");
	}
	next if $used;
	my @nodes = $self->traversal($self->node($node))->get_all;

	my $component   = $self->subgraph(@nodes);
	my @nodes2       = $component->nodes;
	@future{@nodes2} = (1)x@nodes2;
	push(@$components,$component);
    }
    wantarray? @$components: $components;
}

sub by_num {$a <=> $b}

1;

__END__

