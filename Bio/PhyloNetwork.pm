#
# Module for Bio::PhyloNetwork
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Gabriel Cardona <gabriel(dot)cardona(at)uib(dot)es>
#
# Copyright Gabriel Cardona, Gabriel Valiente
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PhyloNetwork - Module to compute with Phylogenetic Networks

=head1 SYNOPSIS

 use Bio::PhyloNetwork;

 # Create a PhyloNetwork object from a eNewick string
 my $net1=Bio::PhyloNetwork->new(
   -eNewick=>'t0:((H1,(H2,l2)),H2); H1:((H3,l1)); H2:((H3,(l3,H1))); H3:(l4);'
 );

 # Print all available data
 print $net1;

 # Rebuild $net1 from its mu_data
 my %mudata=$net1->mudata();
 my $net2=Bio::PhyloNetwork->new(-mudata=>\%mudata,-numleaves=>4);
 print $net2;
 print "d=".$net1->mu_distance($net2)."\n";

 # Get another one and compute distance
 my $net3=Bio::PhyloNetwork->new(
   -eNewick=>'(l2,((l1,(H1,l4)),H1))r; (l3)H1;'
 );
 print "d=".$net1->mu_distance($net3)."\n";

 # ...and find an optimal alignment w.r.t. the Manhattan distance (default)
 my ($weight,%alignment)=$net1->optimal_alignment($net3);
 print "weight:$weight\n";
 foreach my $node1 (keys %alignment) {
   print "$node1 => ".$alignment{$node1}."\n";
 }
 # ...or the Hamming distance

 my ($weightH,%alignmentH)=$net1->optimal_alignment($net3,-metric=>'Hamming');
 print "weight:$weightH\n";
 foreach my $node1 (keys %alignmentH) {
   print "$node1 => ".$alignmentH{$node1}."\n";
 }

 # Test for time consistency of $net1
 if ($net1->is_time_consistent) {
   print "net1 is time consistent\n"
 }
 else {
   print "net1 is not time consistent\n"
 }

 # create a network from the list of edges
 my $net4=Bio::PhyloNetwork->new(-edges=>
   [qw(r s r t s u s c t c t v u b u l3 u b v b v l4 b l2 c l1)]);

 # Test for time consistency of $net3
 if ($net4->is_time_consistent) {
   print "net4 is time consistent\n"
 }
 else {
   print "net4 is not time consistent\n"
 }

 # And print all information on net4
 print $net4;

 # Compute some tripartitions
 my %triparts=$net1->tripartitions();

 # Now these are stored
 print $net1;

 # And can compute the tripartition error
 print "dtr=".$net1->tripartition_error($net3)."\n";

=head1 DESCRIPTION

=head2 Phylogenetic Networks

This is a module to work with phylogenetic networks. Phylogenetic networks
have been studied over the last years as a richer model of the evolutionary
history of sets of organisms than phylogenetic trees, because they take not
only mutation events but also recombination and horizontal gene transfer
events into account.

The natural model for describing the evolutionary
history of a set of sequences under recombination events is a DAG, hence
this package relies on the package Graph::Directed to represent the
underlying graph of a phylogenetic network. We refer the reader to [CRV1,CRV2]
for formal definitions related to phylogenetic networks.

=head2 eNewick description

With this package, phylogenetic networks can be given by its eNewick
string. This description appeared in other packages related to
phylogenetic networks (see [PhyloNet] and [NetGen]); in fact, these two
packages use different descriptions. The Bio::PhyloNetwork
package allows both of them, but uses the second one in its output.

The first approach [PhyloNet] goes as follows: For each hybrid node H, say with
parents u_1,u_2,...,u_k and children v_1,v_2,...v_l: split H in k+1 different
nodes; let each of the first k copies be a child of one of the u_1,...,u_k
(one for each) and have no children (hence we will have k extra leaves);
as for the last copy, let it have no parents and have v_1,...,v_l be its
children. This way we get a forest; each of the trees will be rooted at either
one root of the phylogenetic network or a hybrid node of it; the set of leaves
(of the whole forest) will be the set of leaves of the original network
together with the set of hybrid nodes (each of them repeated as many times
as its in-degree). Then, the eNewick representation of the phylogenetic network
will be the Newick representation of all the trees in the obtained forest,
each of them with its root labeled.

The second approach [NetGen] goes as follows: For each hybrid node H, say with
parents u_1,u_2,...,u_k and children v_1,v_2,...v_l: split H in k different
nodes; let the first copy be a child of u_1 and have all v_1,v_2,...v_l as
its children; let the other copies be child of u_2,...,u_k (one for each)
and have no children. This way, we get a tree whose set of leaves is the
set of leaves of the original network together with the set of hybrid nodes
(possibly repeated). Then the Newick string of the obtained tree (note that
some internal nodes will be labeled and some leaves will be repeated) is
the eNewick string of the phylogenetic network.

For example, consider the network depicted below:

       r
      / \
     /   \
    U     V
   / \   / \
  1   \ /   3
       H
       |
       2

If the first approach is taken, we get the forest:

       r
      / \
     /   \
    U     V
   / \   / \
  1   H H   3
       |
       H
       |
       2

Hence, the eNewick string is '((1,H),(H,3))r; (2)H;'.

As for the second one, one gets the tree:

       r
      / \
     /   \
    U     V
   / \   / \
  1   H |   3
        H
        |
        2

Hence, the eNewick string is '((1,H),((2)H,3))r;'.

Note: when rooting a tree, this package allows the notations
'(subtree,subtree,...)root' as well as 'root:(subtree,subtree,...)', but
the first one is used when writing eNewick strings.

=head2 Tree-child phylogenetic networks

Tree-child (TC) phylogenetic networks are a special class of phylogenetic
networks for which a distance, called mu-distance, is defined [CRV2]
based on certain data (mu-data) associated to every node.
Moreover, this distance extends the
Robinson-Foulds on phylogenetic trees. This package allows testing for a
phylogenetic network if it is TC and computes mu-distances between networks
over the same set of leaves.

Moreover, the mu-data allows one to define the optimal
(in some precise sense) alignment between networks
over the same set of leaves. This package also computes this optimal alignment.

=head2 Tripartitions

Although tripartitions (see [CRV1] and the references therein) do not allow
to define distances, this package outputs tripartitions and computes a weak
form of the tripartition error.

=head2 Time-consistency

Another useful property of Phylogenetic Networks that appears in the literature
is that of time-consistency or real-time hybrids [BSS]. Roughly speaking, a
network admits a temporal representation if it can be drawn in such a way
that tree arcs (those whose end is a tree node) are inclined downwards, while
hybridization arcs (those whose end is a hybrid node) are horizontal.
This package checks for time-consistency and, if so, a temporal representation
is provided.

=head1 AUTHOR

 Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es
 Gabriel Valiente, valiente(at)lsi(dot)upc(dot)edu

=head1 SEE ALSO

=over

=item [CRV1]

G. Cardona, F. Rossello, G. Valiente. Tripartitions do not always
discriminate phylogenetic networks. arXiv:0707.2376v1 [q-bio.PE]

=item [CRV2]

G. Cardona, F. Rossello, G. Valiente. A Distance Measure for
Tree-Child Phylogenetic Networks. Preprint.

=item [NetGen]

M.M. Morin, and B.M.E. Moret. NetGen: generating phylogenetic networks
with diploid hybrids. Bioinformatics 22 (2006), 1921-1923

=item [PhyloNet]

PhyloNet: "Phylogenetic Networks Toolkit".
http://bioinfo.cs.rice.edu/phylonet

=item [BSS]

M. Baroni, C. Semple, and M. Steel. Hybrids in Real
Time. Syst. Biol. 55(1):46-56, 2006

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork;

use strict;
use warnings;

use base qw(Bio::Root::Root);

use Bio::PhyloNetwork::muVector;
use Graph::Directed;
use Bio::TreeIO;
use Bio::Tree::Node;
use IO::String;
use Array::Compare;
use Algorithm::Munkres;

# Creator

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PhyloNetwork();
 Function: Creates a new Bio::PhyloNetwork object
 Returns : Bio::PhyloNetwork
 Args    : none
            OR
           -eNewick => string
            OR
           -graph => Graph::Directed object
            OR
           -edges => reference to an array
            OR
           -tree => Bio::Tree::Tree object
            OR
           -mudata => reference to a hash,
           -leaves => reference to an array
            OR
           -mudata => reference to a hash,
           -numleaves => integer

Returns a Bio::PhyloNetwork object, created according to the data given:

=over 3

=item new()

creates an empty network.

=item new(-eNewick =E<gt> $str)

creates the network whose
Extended Newick representation (see description above) is the string $str.

=item new(-graph =E<gt> $graph)

creates the network with underlying
graph given by the Graph::Directed object $graph

=item new(-tree =E<gt> $tree)

creates a network as a copy of the
Bio::Tree::Tree object in $tree

=item new(-mudata =E<gt> \%mudata, -leaves =E<gt> \@leaves)

creates the network by reconstructing it from its mu-data stored in
\%mudata and with set of leaves in \@leaves.

=item new(-mudata =E<gt> \%mudata, -numleaves =E<gt> $numleaves)

creates the network by reconstructing it from its mu-data stored in
\%mudata and with set of leaves in ("l1".."l$numleaves").

=back

=cut

sub new {
  my ($pkg,@args)=@_;
  my $self=$pkg->SUPER::new(@args);
  my ($eNewick,$edgesR,$leavesR,$numleaves,$graph,$tree,$mudataR)=
    $self->_rearrange([qw(ENEWICK
			  EDGES
			  LEAVES
			  NUMLEAVES
			  GRAPH
			  TREE
			  MUDATA)],@args);
  bless($self,$pkg);

  $self->build_from_eNewick($eNewick) if defined $eNewick;
  $self->build_from_edges(@$edgesR) if defined $edgesR;
  $self->build_from_graph($graph) if defined $graph;
  $self->build_from_tree($tree) if defined $tree;
  if ((! defined $leavesR) && (defined $numleaves)) {
    my @leaves=map {"l$_"} (1..$numleaves);
    $leavesR=\@leaves;
  }
  $self->build_from_mudata($mudataR,$leavesR)
    if ((defined $mudataR) && (defined $leavesR));
  return $self;
}

# Builders

sub build_from_edges {
  my ($self,@edges)=@_;
  my $graph=Graph::Directed->new();
  $graph->add_edges(@edges);
  $self->{graph}=$graph;
  $self->recompute();
  my $labels={};
  foreach my $node ($self->nodes()) {
    $labels->{$node}=$node;
  }
  $self->{labels}=$labels;
}

sub build_from_graph {
  my ($self,$graph)=@_;
  my $graphcp=$graph->copy();
  $self->{graph}=$graphcp;
  $self->recompute();
  my $labels={};
  foreach my $node ($self->nodes()) {
    $labels->{$node}=$node;
  }
  $self->{labels}=$labels;
}

my $_eN_index;

sub build_from_eNewick {
  my ($self,$string)=@_;
  $_eN_index=0;
  my $graph=Graph::Directed->new();
  my $labels={};
  my @blocks=split(/; */,$string);
  foreach my $block (@blocks) {
    my ($rt,$str)=get_root_and_subtree($block);
    my ($rtlbl,$rttype,$rtid,$rtlng)=get_label_type_id_length($rt);
    process_block($graph,$labels,$block,$rtid);
    $labels->{$rtid}=$rtlbl.'';
  }
  $self->{graph}=$graph;
  $self->{labels}=$labels;
  $self->recompute();
}

sub process_block {
  my ($graph,$labels,$block,$rtid)=@_;
  my ($rt,$str)=get_root_and_subtree($block);
  my @substrs=my_split($str);
  foreach my $substr (@substrs) {
    my ($subrt,$subblock)=get_root_and_subtree($substr);
    my ($subrtlbl,$subrttype,$subrtid,$subrtlng)=
      get_label_type_id_length($subrt);
    if (! $subrtlng eq '') {
      $graph->add_weighted_edges($rtid,$subrtid,$subrtlng);
    }
    else {
      $graph->add_edges($rtid,$subrtid);
    }
    if (! $subrttype eq '') {
      $graph->set_edge_attribute($rtid,$subrtid,'type',$subrttype);
    }
    $subrtlbl.='';
#    if (! $subrtlbl eq '') {
    if ((! defined $labels->{$subrtid})||($labels->{$subrtid} eq '')){
      $labels->{$subrtid}=$subrtlbl;
    } elsif (( $labels->{$subrtid} ne $subrtlbl )&&($subrtlbl ne '')) {
      # error
      die("Different labels for the same node (".$labels->{$subrtid}." and $subrtlbl)");
    }
#    }
    if ($subblock ne "") {
      process_block($graph,$labels,$subblock,$subrtid);
    }
  }
}

sub get_root_and_subtree {
  my ($block)=@_;
  my ($rt,$str)=("","");
#  ($rt,$str)=split(/:|=/,$block);
  ($rt,$str)=split(/=/,$block);
  if ($rt eq $block) {
    # try to look for root label at the end
    my $pos=length($rt)-1;
    while ((substr($rt,$pos,1) ne ")") && ($pos >=0)) {
      $pos--;
    }
    $rt=substr($block,$pos+1,length($block)-$pos);
    $str=substr($block,0,$pos+1);
  }
  $rt=trim($rt);
  $str=trim($str);
  return ($rt,$str);
}

sub get_label_type_id_length {
  my ($string) = @_;
  $string.='';
#  print "$string\n";
  if (index($string,'#')==-1) {
    # no hybrid
    my ($label,$length)=split(':',$string);
    $label.='';
    my $id;
    if ((! defined $label) || ($label eq '')) {
      # create id
      $_eN_index++;
      $id="T$_eN_index";
    } else {
      $id=$label;
    }
    return ($label,'',$id,$length);
  }
  else {
    # hybrid
    my ($label,$string2)=split('#',$string);
    my ($typeid,$length)=split(':',$string2);
    my $type=$typeid;
    $type =~ s/\d//g;
    my $id=$typeid;
    $id =~ s/\D//g;
    return ($label,$type,'#'.$id,$length);
  }
}

sub trim
{
  my ($string) = @_;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

sub my_split {
  my ( $string ) = @_;
  my $temp="";
  my @substrings;
  my $level=1;
  for my $i ( 1 .. length( $string ) ) {
    my $char=substr($string,$i,1);
    if ($char eq "(") {
      $level++;
    }
    if ($char eq ")") {
      if ($level==1) {
      	push @substrings, $temp;
	$temp="";
      }
      $level--;
    }
    if (($char eq ",") && ($level==1)) {
      	push @substrings, $temp;
	$temp="";
	$char="";
    }
    $temp = $temp.$char;
  }
  return @substrings;
}

sub build_from_mudata {
  my ($self,$mus,$leavesR)=@_;
  my $graph=Graph::Directed->new();
  my @nodes=keys %{$mus};
  my @leaves=@{$leavesR};

  my %seen;
  my @internal;

  @seen{@leaves} = ();

  foreach my $node (@nodes) {
    push(@internal, $node) unless exists $seen{$node};
  }

  @internal=sort {$mus->{$b} <=> $mus->{$a} } @internal;
  @nodes=(@internal,@leaves);
  my $numnodes=@nodes;
  for (my $i=0;$i<$numnodes;$i++) {
    my $mu=$mus->{$nodes[$i]};
    my $j=$i+1;
    while ($mu->is_positive() && $j<$numnodes) {
      if ($mu->geq_poset($mus->{$nodes[$j]})) {
	$graph->add_edges(($nodes[$i],$nodes[$j]));
	$mu = $mu - $mus->{$nodes[$j]};
      }
      $j++;
    }
  }
  $self->build_from_graph($graph);
}

# sub relabel_tree {
#   my ($tree)=@_;
#   my $i=1;
#   my $j=1;
#   my $root=$tree->get_root_node();
#   foreach my $node ($tree->get_nodes()) {
#     if ($node == $root) {
#       $node->{'_id'}="r";
#     }
#     elsif (! $node->is_Leaf) {
#       $node->{'_id'}="t$i";
#       $i++;
#     }
#     else {
#       if ($node->{'_id'} eq "") {
# 	$node->{'_id'}="l$j";
# 	$j++;
#       }
#     }
#   }
#   return $tree;
# }

# sub build_subtree {
#   my ($graph,$root)=@_;
#   foreach my $child ($root->each_Descendent) {
#     $graph->add_edge($root->id,$child->id);
#     $graph=build_subtree($graph,$child);
#   }
#   return $graph;
# }

sub build_from_tree {
  my ($self,$tree)=@_;
#  relabel_tree($tree);
#  my $treeroot=$tree->get_root_node;
#  my $graph=Graph::Directed->new();
#  $graph=build_subtree($graph,$treeroot);
#  $self->build_from_graph($graph);
  my $str;
  my $io=IO::String->new($str);
  my $treeio=Bio::TreeIO->new(-format => 'newick', -fh => $io);
  $treeio->write_tree($tree);
#  print "intern: $str\n";
  $self->build_from_eNewick($str);
}

sub recompute {
  my ($self)=@_;
  $self->throw("Graph is not DAG:".$self->{graph})
    unless $self->{graph}->is_dag();
  my @leaves=$self->{graph}->successorless_vertices();
  @leaves=sort @leaves;
  my $numleaves=@leaves;
  my @roots=$self->{graph}->predecessorless_vertices();
  my $numroots=@roots;
  #$self->throw("Graph is not rooted") unless ($numroots == 1);
  my @nodes=$self->{graph}->vertices();
  @nodes=sort @nodes;
  my $numnodes=@nodes;
  foreach my $node (@nodes) {
    if (! defined $self->{labels}->{$node}) {
      $self->{labels}->{$node}='';
    }
  }
  $self->{leaves}=\@leaves;
  $self->{numleaves}=$numleaves;
  $self->{roots}=\@roots;
  $self->{numroots}=$numroots;
  $self->{nodes}=\@nodes;
  $self->{numnodes}=$numnodes;
  $self->{mudata}={};
  $self->{h}={};
  $self->compute_height();
  $self->compute_mu();
  return $self;
}

# Hybridizing

sub is_attackable {
  my ($self,$u1,$v1,$u2,$v2)=@_;
  if ( $self->is_hybrid_node($v1) ||
       $self->is_hybrid_node($v2) ||
       $self->graph->is_reachable($v2,$u1) ||
       (($u1 eq $u2)&&($v1 eq $v2)) ||
       (! scalar grep {($_ ne $v2) && ($self->is_tree_node($_))}
	$self->graph->successors($u2)))
    {
      return 0;
    }
  return 1;
}

sub do_attack {
  my ($self,$u1,$v1,$u2,$v2,$lbl)=@_;
  my $graph=$self->{graph};
  $graph->delete_edge($u1,$v1);
  $graph->delete_edge($u2,$v2);
  $graph->add_edge($u1,"T$lbl");
  $graph->add_edge("T$lbl",$v1);
  $graph->add_edge($u2,"#H$lbl");
  $graph->add_edge("#H$lbl",$v2);
  $graph->add_edge("T$lbl","#H$lbl");
  $self->build_from_graph($graph);
}


# Computation of mu-data

sub compute_mu {
  my ($self)=@_;
  my $graph=$self->{graph};
  my $mudata=$self->{mudata};
  my @leaves=@{$self->{leaves}};
  my $numleaves=$self->{numleaves};
  for (my $i=0;$i<$numleaves;$i++) {
    my $vec=Bio::PhyloNetwork::muVector->new($numleaves);
    $vec->[$i]=1;
    $mudata->{$leaves[$i]}=$vec;
  }
  my $h=1;
  while (my @nodes=grep {$self->{h}->{$_} == $h} @{$self->{nodes}} )
    {
      foreach my $u (@nodes) {
	my $vec=Bio::PhyloNetwork::muVector->new($numleaves);
	foreach my $son ($graph->successors($u)) {
	  $vec+=$mudata->{$son};
	}
	$mudata->{$u}=$vec;
      }
      $h++;
    }
}

sub compute_height {
  my ($self)=@_;
  my $graph=$self->{graph};
  my @leaves=@{$self->{leaves}};
  foreach my $leaf (@leaves) {
    $self->{h}->{$leaf}=0;
  }
  my $h=0;
  while (my @nodes=grep {(defined $self->{h}->{$_})&&($self->{h}->{$_} == $h)}
	 @{$self->{nodes}} )
    {
    foreach my $node (@nodes) {
      foreach my $parent ($graph->predecessors($node)) {
	$self->{h}->{$parent}=$h+1;
      }
    }
    $h++;
  }
}

# Tests

=head2 is_leaf

 Title   : is_leaf
 Usage   : my $b=$net->is_leaf($u)
 Function: tests if $u is a leaf in $net
 Returns : boolean
 Args    : scalar

=cut

sub is_leaf {
  my ($self,$node)=@_;
  if ($self->{graph}->out_degree($node) == 0) {return 1;}
  return 0;
}

=head2 is_root

 Title   : is_root
 Usage   : my $b=$net->is_root($u)
 Function: tests if $u is the root of $net
 Returns : boolean
 Args    : scalar

=cut

sub is_root {
  my ($self,$node)=@_;
  if ($self->{graph}->in_degree($node) == 0) {return 1;}
  return 0;
}

=head2 is_tree_node

 Title   : is_tree_node
 Usage   : my $b=$net->is_tree_node($u)
 Function: tests if $u is a tree node in $net
 Returns : boolean
 Args    : scalar

=cut

sub is_tree_node {
  my ($self,$node)=@_;
  if ($self->{graph}->in_degree($node) <= 1) {return 1;}
  return 0;
}

=head2 is_hybrid_node

 Title   : is_hybrid_node
 Usage   : my $b=$net->is_hybrid_node($u)
 Function: tests if $u is a hybrid node in $net
 Returns : boolean
 Args    : scalar

=cut

sub is_hybrid_node {
  my ($self,$node)=@_;
  if ($self->{graph}->in_degree($node) > 1) {return 1;}
  return 0;
}

sub has_tree_child {
  # has_tree_child(g,u) returns 1 if u has a tree child in graph g
  # and 0 otherwise
  my $g=shift(@_);
  my $node=shift(@_);
  my @Sons=$g->successors($node);
  foreach my $son (@Sons) {
    if ($g->in_degree($son)==1) {
      return 1;
    }
  }
  return 0;
}

=head2 is_tree_child

 Title   : is_tree_child
 Usage   : my $b=$net->is_tree_child()
 Function: tests if $net is a Tree-Child phylogenetic network
 Returns : boolean
 Args    : Bio::PhyloNetwork

=cut

sub is_tree_child {
  my ($self)=@_;
  if (defined $self->{is_tree_child}) {
    return $self->{is_tree_child};
  }
  $self->{is_tree_child}=0;
  my $graph=$self->{graph};
  foreach my $node (@{$self->{nodes}}) {
    return 0 unless ($graph->out_degree($node)==0 ||
		     has_tree_child($graph,$node));
  }
  $self->{is_tree_child}=1;
  return 1;
}

# Accessors

=head2 nodes

 Title   : nodes
 Usage   : my @nodes=$net->nodes()
 Function: returns the set of nodes of $net
 Returns : array
 Args    : none

=cut

sub nodes {
  my ($self)=@_;
  return @{$self->{nodes}};
}

=head2 leaves

 Title   : leaves
 Usage   : my @leaves=$net->leaves()
 Function: returns the set of leaves of $net
 Returns : array
 Args    : none

=cut

sub leaves {
  my ($self)=@_;
  return @{$self->{leaves}};
}

=head2 roots

 Title   : roots
 Usage   : my @roots=$net->roots()
 Function: returns the set of roots of $net
 Returns : array
 Args    : none

=cut

sub roots {
  my ($self)=@_;
  return @{$self->{roots}};
}

=head2 internal_nodes

 Title   : internal_nodes
 Usage   : my @internal_nodes=$net->internal_nodes()
 Function: returns the set of internal nodes of $net
 Returns : array
 Args    : none

=cut

sub internal_nodes {
  my ($self)=@_;
  return grep {! $self->is_leaf($_)} $self->nodes();
}

=head2 tree_nodes

 Title   : tree_nodes
 Usage   : my @tree_nodes=$net->tree_nodes()
 Function: returns the set of tree nodes of $net
 Returns : array
 Args    : none

=cut

sub tree_nodes {
  my ($self)=@_;
  return grep {$self->is_tree_node($_)} $self->nodes();
}

=head2 hybrid_nodes

 Title   : hybrid_nodes
 Usage   : my @hybrid_nodes=$net->hybrid_nodes()
 Function: returns the set of hybrid nodes of $net
 Returns : array
 Args    : none

=cut

sub hybrid_nodes {
  my ($self)=@_;
  return grep {$self->is_hybrid_node($_)} $self->nodes();
}

=head2 graph

 Title   : graph
 Usage   : my $graph=$net->graph()
 Function: returns the underlying graph of $net
 Returns : Graph::Directed
 Args    : none

=cut

sub graph {
  my ($self)=@_;
  return $self->{graph};
}

=head2 edges

 Title   : edges
 Usage   : my @edges=$net->edges()
 Function: returns the set of edges of $net
 Returns : array
 Args    : none

Each element in the array is an anonimous array whose first element is the
head of the edge and the second one is the tail.

=cut

sub edges {
  my ($self)=@_;
  return $self->{graph}->edges();
}

=head2 tree_edges

 Title   : tree_edges
 Usage   : my @tree_edges=$net->tree_edges()
 Function: returns the set of tree edges of $net
           (those whose tail is a tree node)
 Returns : array
 Args    : none

=cut

sub tree_edges {
  my ($self)=@_;
  return grep {$self->is_tree_node($_->[1])} $self->edges();
}

=head2 hybrid_edges

 Title   : hybrid_edges
 Usage   : my @hybrid_edges=$net->hybrid_edges()
 Function: returns the set of hybrid edges of $net
           (those whose tail is a hybrid node)
 Returns : array
 Args    : none

=cut

sub hybrid_edges {
  my ($self)=@_;
  return grep {$self->is_hybrid_node($_->[1])} $self->edges();
}

=head2 explode

 Title   : explode
 Usage   : my @trees=$net->explode()
 Function: returns the representation of $net by a set of
           Bio::Tree:Tree objects
 Returns : array
 Args    : none

=cut

sub explode {
  my ($self)=@_;
  my @trees;
  $self->explode_rec(\@trees);
  return @trees;
}

sub explode_rec {
  my ($self,$trees)=@_;
  my @h = $self->hybrid_nodes;
  if (scalar @h) {
    my $v = shift @h;
    for my $u ($self->{graph}->predecessors($v)) {
      $self->{graph}->delete_edge($u,$v);
      $self->explode_rec($trees);
      $self->{graph}->add_edge($u,$v);
    }
  } else {
    my $io = IO::String->new($self->eNewick);
    my $treeio = Bio::TreeIO->new(-format => 'newick', -fh => $io);
    my $tree = $treeio->next_tree;
    $tree->contract_linear_paths;
    push @{$trees}, $tree;
  }
}

=head2 mudata

 Title   : mudata
 Usage   : my %mudata=$net->mudata()
 Function: returns the representation of $net by its mu-data
 Returns : hash
 Args    : none

$net-E<gt>mudata() returns a hash with keys the nodes of $net and each value is a
muVector object holding its mu-vector.

=cut

sub mudata {
  my ($self)=@_;
  return %{$self->{mudata}};
}

sub mudata_node {
  my ($self,$u)=@_;
  return $self->{mudata}{$u};
}

=head2 heights

 Title   : heights
 Usage   : my %heights=$net->heights()
 Function: returns the heights of the nodes of $net
 Returns : hash
 Args    : none

$net-E<gt>heights() returns a hash with keys the nodes of $net and each value
is its height.

=cut

sub heights {
  my ($self)=@_;
  return %{$self->{h}};
}

sub height_node {
  my ($self,$u)=@_;
  return $self->{h}{$u};
}

=head2 mu_distance

 Title   : mu_distance
 Usage   : my $dist=$net1->mu_distance($net2)
 Function: Computes the mu-distance between the networks $net1 and $net2 on
           the same set of leaves
 Returns : scalar
 Args    : Bio::PhyloNetwork

=cut

sub mu_distance {
  my ($net1,$net2)=@_;
  my @nodes1=$net1->nodes;
  my @nodes2=$net2->nodes;
  my $comp = Array::Compare->new;
  $net1->throw("Cannot compare phylogenetic networks on different set of leaves")
    unless $comp->compare($net1->{leaves},$net2->{leaves});
  $net1->warn("Not a tree-child phylogenetic network")
    unless $net1->is_tree_child();
  $net2->warn("Not a tree-child phylogenetic network")
    unless $net2->is_tree_child();
  my @leaves=@{$net1->{leaves}};
  my %matched1;
  my %matched2;
  OUTER: foreach my $node1 (@nodes1) {
    foreach my $node2 (@nodes2) {
      if (
	  (! exists $matched1{$node1}) && (! exists $matched2{$node2}) &&
	  ($net1->{mudata}{$node1} == $net2->{mudata}{$node2})
	 ) {
	$matched1{$node1}=$node2;
	$matched2{$node2}=$node1;
	next OUTER;
      }
    }
  }
  return (scalar @nodes1)+(scalar @nodes2)-2*(scalar keys %matched1);
}

=head2 mu_distance_generalized

 Title   : mu_distance_generalized
 Usage   : my $dist=$net1->mu_distance($net2)
 Function: Computes the mu-distance between the topological restrictions of
           networks $net1 and $net2 on its common set of leaves
 Returns : scalar
 Args    : Bio::PhyloNetwork

=cut

sub mu_distance_generalized {
  my ($net1,$net2)=@_;
  my ($netr1,$netr2)=$net1->topological_restriction($net2);
  return $netr1->mu_distance($netr2);
}

# mudata_string (code mu_data in a string; useful for isomorphism testing)

sub mudata_string_node {
  my ($self,$u)=@_;
  return $self->{mudata}->{$u}->display();
}

sub mudata_string {
  my ($self)=@_;
  return $self->{mudata_string} if defined $self->{mudata_string};
  my @internal=$self->internal_nodes;
  my $mus=$self->{mudata};
  @internal=sort {$mus->{$b} <=> $mus->{$a} } @internal;
  my $str="";
  foreach my $node (@internal) {
    $str=$str.$self->mudata_string_node($node);
  }
  $self->{mudata_string}=$str;
  return $str;
}

sub is_mu_isomorphic {
  my ($net1,$net2)=@_;
  return ($net1->mudata_string() eq $net2->mudata_string());
}

# tripartitions

sub compute_tripartition_node {
  my ($self,$u)=@_;
  $self->warn("Cannot compute tripartitions on unrooted networks. Will assume one at random")
    unless ($self->{numroots} == 1);
  my $root=$self->{roots}->[0];
  my $graph=$self->{graph};
  my $graphPruned=$graph->copy();
  $graphPruned->delete_vertex($u);
  my $tripartition="";
  foreach my $leaf (@{$self->{leaves}}) {
    my $type;
    if ($graph->is_reachable($u,$leaf)) {
      if ($graphPruned->is_reachable($root,$leaf)) {$type="B";}
      else {$type="A";}
    }
    else {$type="C";}
    $tripartition .= $type;
  }
  $self->{tripartitions}->{$u}=$tripartition;
}

sub compute_tripartitions {
  my ($self)=@_;
  foreach my $node (@{$self->{nodes}}) {
    $self->compute_tripartition_node($node);
  }
}

=head2 tripartitions

 Title   : tripartitions
 Usage   : my %tripartitions=$net->tripartitions()
 Function: returns the set of tripartitions of $net
 Returns : hash
 Args    : none

$net-E<gt>tripartitions() returns a hash with keys the nodes of $net and each value
is a string representing the tripartition of the leaves induced by the node.
A string "BCA..." associated with a node u (e.g.) means, the first leaf is in
the set B(u), the second one in C(u), the third one in A(u), and so on.

=cut

sub tripartitions {
  my ($self)=@_;
  $self->compute_tripartitions() unless defined $self->{tripartitions};
  return %{$self->{tripartitions}};
}

# to do: change to tri_distance and test for TC and time-cons

sub tripartition_error {
  my ($net1,$net2)=@_;
  my $comp = Array::Compare->new;
  $net1->throw("Cannot compare phylogenetic networks on different set of leaves")
    unless $comp->compare($net1->{leaves},$net2->{leaves});
  $net1->warn("Not a tree-child phylogenetic network")
    unless $net1->is_tree_child();
  $net2->warn("Not a tree-child phylogenetic network")
    unless $net2->is_tree_child();
  $net1->warn("Not a time-consistent network")
    unless $net1->is_time_consistent();
  $net2->warn("Not a time-consistent network")
    unless $net2->is_time_consistent();
  $net1->compute_tripartitions() unless defined $net1->{tripartitions};
  $net2->compute_tripartitions() unless defined $net2->{tripartitions};
  my @edges1=$net1->{graph}->edges();
  my @edges2=$net2->{graph}->edges();
  my ($FN,$FP)=(0,0);
  foreach my $edge1 (@edges1) {
    my $matched=0;
    foreach my $edge2 (@edges2) {
      if ($net1->{tripartitions}->{$edge1->[1]} eq
	  $net2->{tripartitions}->{$edge2->[1]}) {
	$matched=1;
	last;
      }
    }
    if (! $matched) {$FN++;}
  }
  foreach my $edge2 (@edges2) {
    my $matched=0;
    foreach my $edge1 (@edges1) {
      if ($net1->{tripartitions}->{$edge1->[1]} eq
	  $net2->{tripartitions}->{$edge2->[1]}) {
	$matched=1;
	last;
      }
    }
    if (! $matched) {$FP++;}
  }
  return ($FN/(scalar @edges1)+$FP/(scalar @edges2))/2;
}

# Time-consistency

# to do: add weak time consistency

=head2 is_time_consistent

 Title   : is_time_consistent
 Usage   : my $b=$net->is_time_consistent()
 Function: tests if $net is (strong) time-consistent
 Returns : boolean
 Args    : none

=cut

sub is_time_consistent {
  my ($self)=@_;
  $self->compute_temporal_representation()
    unless exists $self->{has_temporal_representation};
  return $self->{has_temporal_representation};
}

=head2 temporal_representation

 Title   : temporal_representation
 Usage   : my %time=$net->temporal_representation()
 Function: returns a hash containing a temporal representation of $net, or 0
           if $net is not time-consistent
 Returns : hash
 Args    : none

=cut

sub temporal_representation {
  my ($self)=@_;
  if ($self->is_time_consistent) {
    return %{$self->{temporal_representation}};
  }
  return 0;
}

sub compute_temporal_representation {
  my ($self)=@_;
  my $quotient=Graph::Directed->new();
  my $classes=find_classes($self);
  my %repr;
  map {$repr{$_}=$classes->{$_}[0]} $self->nodes();
  foreach my $e ($self->tree_edges()) {
    $quotient->add_edge($repr{$e->[0]},$repr{$e->[1]});
  }
  my %temp;
  my $depth=0;
  while ($quotient->vertices()) {
    if (my @svs=$quotient->predecessorless_vertices()) {
      foreach my $sv (@svs) {
	$temp{$sv}=$depth;
      }
      $quotient->delete_vertices(@svs);
    } else {
      return 0;
    }
    $depth++;
  }
  foreach my $node (@{$self->{nodes}}) {
    $temp{$node}=$temp{$repr{$node}}
  }
  $self->{temporal_representation}=\%temp;
  $self->{has_temporal_representation}=1;
}

sub find_classes {
  my ($self)=@_;
  my $classes={};
  map {$classes->{$_}=[$_]} $self->nodes();
  foreach my $e ($self->hybrid_edges()) {
    $classes=join_classes($classes,$e->[0],$e->[1]);
  }
  return $classes;
}

sub join_classes {
  my ($classes,$u,$v)=@_;
  my @clu=@{$classes->{$u}};
  my @clv=@{$classes->{$v}};
  my @cljoin=(@clu,@clv);
  map {$classes->{$_}=\@cljoin} @cljoin;
  return $classes;
}

# alignment

=head2 contract_elementary


 Title   : contract_elementary
 Usage   : my ($contracted,$blocks)=$net->contract_elementary();
 Function: Returns the network $contracted, obtained by contracting elementary
           paths of $net into edges. The reference $blocks points to a hash
           where, for each node of $contracted, gives the corresponding nodes
           of $net that have been deleted.
 Returns : Bio::PhyloNetwork,reference to hash
 Args    : none

=cut

sub contract_elementary {
  my ($self)=@_;

  my $contracted=$self->graph->copy();
  my @nodes=$self->nodes();
  my $mus=$self->{mudata};
  my $hs=$self->{h};
  my %blocks;
  foreach my $u (@nodes) {
    $blocks{$u}=[$u];
  }
  my @elementary=grep { $contracted->out_degree($_) == 1} $self->tree_nodes();
  @elementary=sort {$mus->{$b} <=> $mus->{$a} ||
			 $hs->{$b} <=> $hs->{$a}} @elementary;
  foreach my $elem (@elementary) {
    my @children=$contracted->successors($elem);
    my $child=$children[0];
    if ($contracted->in_degree($elem) == 1) {
      my @parents=$contracted->predecessors($elem);
      my $parent=$parents[0];
      $contracted->add_edge($parent,$child);
    }
    $contracted->delete_vertex($elem);
    my @blch=@{$blocks{$child}};
    my @blem=@{$blocks{$elem}};
    $blocks{$child}=[@blem,@blch];
    delete $blocks{$elem};
  }
  my $contr=Bio::PhyloNetwork->new(-graph => $contracted);
  return $contr,\%blocks;
}

=head2 optimal_alignment

 Title   : optimal_alignment
 Usage   : my ($weight,$alignment,$wgts)=$net->optimal_alignment($net2)
 Function: returns the total weight of an optimal alignment,
           the alignment itself, and partial weights
           between the networks $net1 and $net2 on the same set of leaves.
           An optional argument allows one to use the Manhattan (default) or the
           Hamming distance between mu-vectors.
 Returns : scalar,reference to hash,reference to hash
 Args    : Bio::PhyloNetwork,
           -metric => string (optional)

Supported strings for the -metric parameter are 'Manhattan' or 'Hamming'.

=cut

sub optimal_alignment {
  my ($net1,$net2,%params)=@_;

  my ($net1cont,$blocks1)=contract_elementary($net1);
  my ($net2cont,$blocks2)=contract_elementary($net2);
  my ($wc,$alignc,$weightc)=
    optimal_alignment_noelementary($net1cont,$net2cont,%params);
  my %alignment=();
  my $totalweigth=0;
  my %weigths=();
  foreach my $u1 (keys %$alignc) {
    my $u2=$alignc->{$u1};
    my @block1=@{$blocks1->{$u1}};
    my @block2=@{$blocks2->{$u2}};
    while (@block1 && @block2) {
      my $u1dc=pop @block1;
      my $u2dc=pop @block2;
      $alignment{$u1dc}=$u2dc;
      $weigths{$u1dc}=$weightc->{$u1};
      $totalweigth+=$weigths{$u1dc};
    }
  }
  return $totalweigth,\%alignment,\%weigths;
}

sub optimal_alignment_noelementary {
  my ($net1,$net2,%params)=@_;

  my $comp = Array::Compare->new;
  $net1->throw("Cannot align phylogenetic networks on different set of leaves")
    unless $comp->compare($net1->{leaves},$net2->{leaves});
  my $distance;
  if ((defined $params{-metric})and ($params{-metric} eq 'Hamming')) {
    $distance='Hamming';
  } else {
    $distance='Manhattan';
  }
  my $numleaves=$net1->{numleaves};
  my @nodes1=$net1->internal_nodes();
  my @nodes2=$net2->internal_nodes();
  my $numnodes1=@nodes1;
  my $numnodes2=@nodes2;
  my @matrix=();
  for (my $i=0;$i<$numnodes1;$i++) {
    my @row=();
    for (my $j=0;$j<$numnodes2;$j++) {
      push @row,weight($net1,$nodes1[$i],$net2,$nodes2[$j],$distance);
    }
    push @matrix,\@row;
  }
  my @alignment=();
  Algorithm::Munkres::assign(\@matrix,\@alignment);
  my %alignmenthash;
  my %weighthash;
  my $totalw=0;
  foreach my $leaf (@{$net1->{leaves}}) {
    $alignmenthash{$leaf}=$leaf;
    $weighthash{$leaf}=0;
  }
  for (my $i=0;$i<$numnodes1;$i++) {
    if (defined $nodes2[$alignment[$i]]) {
      $alignmenthash{$nodes1[$i]}=$nodes2[$alignment[$i]];
      $weighthash{$nodes1[$i]}=$matrix[$i][$alignment[$i]];
      $totalw += $matrix[$i][$alignment[$i]];
    }
  }
  return $totalw,\%alignmenthash,\%weighthash;
 }

=head2 optimal_alignment_generalized

 Title   : optimal_alignment_generalized
 Usage   : my ($weight,%alignment)=$net->optimal_alignment_generalized($net2)
 Function: returns the wieght of an optimal alignment, and the alignment itself,
           between the topological restriction of the networks $net1 and $net2
           on the set of common leaves.
           An optional argument allows one to use the Manhattan (default) or the
           Hamming distance between mu-vectors.
 Returns : scalar,hash
 Args    : Bio::PhyloNetwork,
           -metric => string (optional)

Supported strings for the -metric parameter are 'Manhattan' or 'Hamming'.

=cut

sub optimal_alignment_generalized {
  my ($net1,$net2,%params)=@_;
  my ($netr1,$netr2)=$net1->topological_restriction($net2);
  return $netr1->optimal_alignment($netr2,%params);
}

sub weight {
  my ($net1,$v1,$net2,$v2,$distance)=@_;
  my $w;
  if (! defined $distance) {
    $distance='Manhattan';
  }
  if ($distance eq 'Hamming') {
    $w=$net1->{mudata}->{$v1}->hamming($net2->{mudata}->{$v2});
  } else {
    $w=$net1->{mudata}->{$v1}->manhattan($net2->{mudata}->{$v2});
  }
  if (($net1->is_tree_node($v1) && $net2->is_hybrid_node($v2)) ||
      ($net2->is_tree_node($v2) && $net1->is_hybrid_node($v1))
     )
    {
      $w +=1/(2*$net1->{numleaves});
    }
  return $w;
}


=head2 topological_restriction

 Title   : topological_restriction
 Usage   : my ($netr1,$netr2)=$net1->topological_restriction($net2)
 Function: returns the topological restriction of $net1 and $net2 on its
           common set of leaves
 Returns : Bio::PhyloNetwork, Bio::PhyloNetwork
 Args    : Bio::PhyloNetwork

=cut

sub topological_restriction {
  my ($net1,$net2)=@_;

  my @leaves1=$net1->leaves();
  my @leaves2=$net2->leaves();
  my $numleaves1=scalar @leaves1;
  my $numleaves2=scalar @leaves2;
  my %position1;
  for (my $i=0; $i<$numleaves1; $i++) {
    $position1{$leaves1[$i]}=$i;
  }
  my %position2;
  my @commonleaves=();
  for (my $j=0; $j<$numleaves2; $j++) {
    if (defined $position1{$leaves2[$j]}) {
      push @commonleaves,$leaves2[$j];
      $position2{$leaves2[$j]}=$j;
    }
  }
  my $graphred1=$net1->{graph}->copy();
  my $graphred2=$net2->{graph}->copy();
 OUTER1:
  foreach my $u ($graphred1->vertices()) {
    my $mu=$net1->mudata_node($u);
    foreach my $leaf (@commonleaves) {
      if ($mu->[$position1{$leaf}]>0) {
	next OUTER1;
      }
    }
    $graphred1->delete_vertex($u);
  }
 OUTER2:
  foreach my $u ($graphred2->vertices()) {
    my $mu=$net2->mudata_node($u);
    foreach my $leaf (@commonleaves) {
      if ($mu->[$position2{$leaf}]>0) {
	next OUTER2;
      }
    }
    $graphred2->delete_vertex($u);
  }
  my $netr1=Bio::PhyloNetwork->new(-graph => $graphred1);
  my $netr2=Bio::PhyloNetwork->new(-graph => $graphred2);
  return ($netr1,$netr2);
}

# Functions for eNewick representation

=head2 eNewick

 Title   : eNewick
 Usage   : my $str=$net->eNewick()
 Function: returns the eNewick representation of $net without labeling
           internal tree nodes
 Returns : string
 Args    : none

=cut

sub eNewick {
  my ($self)=@_;
  my $str="";
  my $seen={};
  foreach my $root ($self->roots()) {
    $str=$str.$self->eNewick_aux($root,$seen,undef)."; ";
  }
  return $str;
}

sub eNewick_aux {
  my ($self,$node,$seen,$parent)=@_;
  my $str='';
  if ($self->is_leaf($node) ||
      (defined $seen->{$node}) )
    {
      $str=make_label($self,$parent,$node);
    }
  else {
    $seen->{$node}=1;
    my @sons=$self->{graph}->successors($node);
    $str="(";
    foreach my $son (@sons) {
      $str=$str.$self->eNewick_aux($son,$seen,$node).",";
    }
    chop($str);
    $str.=")".make_label($self,$parent,$node);
  }
  return $str;
}

sub make_label {
  my ($self,$parent,$node)=@_;
  my $str='';
  if ($self->is_hybrid_node($node)) {
    my $lbl=$self->{labels}->{$node};
    if ($lbl =~ /#/) {
      $lbl='';
    }
    $str.=$lbl; #$self->{labels}->{$node};
    $str.='#';
    if ((defined $parent) &&
	($self->graph->has_edge_attribute($parent,$node,'type'))) {
      $str.=$self->graph->get_edge_attribute($parent,$node,'type');
    }
    $str.=substr $node,1;
  } else {
    $str.=$self->{labels}->{$node};
  }
  if ((defined $parent) &&
      ($self->graph->has_edge_weight($parent,$node))) {
    $str.=":".$self->graph->get_edge_weight($parent,$node);
  }
  return $str;
}

=head2 eNewick_full

 Title   : eNewick_full
 Usage   : my $str=$net->eNewick_full()
 Function: returns the eNewick representation of $net labeling
           internal tree nodes
 Returns : string
 Args    : none

=cut

sub eNewick_full {
  my ($self)=@_;
  my $str="";
  my $seen={};
  foreach my $root ($self->roots()) {
    $str=$str.$self->eNewick_full_aux($root,$seen,undef)."; ";
  }
  return $str;
}

sub eNewick_full_aux {
  my ($self,$node,$seen,$parent)=@_;
  my $str='';
  if ($self->is_leaf($node) ||
      (defined $seen->{$node}) )
    {
      $str=make_label_full($self,$parent,$node);
    }
  else {
    $seen->{$node}=1;
    my @sons=$self->{graph}->successors($node);
    $str="(";
    foreach my $son (@sons) {
      $str=$str.$self->eNewick_full_aux($son,$seen,$node).",";
    }
    chop($str);
    $str.=")".make_label_full($self,$parent,$node);
  }
  return $str;
}

sub make_label_full {
  my ($self,$parent,$node)=@_;
  my $str='';
  if ($self->is_hybrid_node($node)) {
    my $lbl=$self->{labels}->{$node};
    if ($lbl =~ /#/) {
      $lbl='';
    }
    $str.=$lbl; #$self->{labels}->{$node};
    $str.='#';
    if ((defined $parent) &&
	($self->graph->has_edge_attribute($parent,$node,'type'))) {
      $str.=$self->graph->get_edge_attribute($parent,$node,'type');
    }
    $str.=substr $node,1;
  } else {
    if ((defined $self->{labels}->{$node})&&($self->{labels}->{$node} ne '')) {
      $str.=$self->{labels}->{$node};
    }
    else {
      $str.=$node;
    }
  }
  if ((defined $parent) &&
      ($self->graph->has_edge_weight($parent,$node))) {
    $str.=":".$self->graph->get_edge_weight($parent,$node);
  }
  return $str;
}

# sub eNewick_full {
#   my ($self)=@_;
#   my $str="";
#   my $seen={};
#   foreach my $root ($self->roots()) {
#     $str=$str.$self->eNewick_full_aux($root,$seen,undef)."; ";
#   }
#   return $str;
# }

# sub eNewick_full_aux {
#   my ($self,$node,$seen,$parent)=@_;
#   my $str;
#   if ($self->is_leaf($node) ||
#       (defined $seen->{$node}) )
#     {
#       if ($self->is_hybrid_node($node)) {
# 	my $tag=substr $node,1;
# 	if ((defined $parent) &&
# 	    ($self->graph->has_edge_attribute($parent,$node,'type'))) {
# 	  $str='#'.$self->graph->get_edge_attribute($parent,$node,'type').$tag;
# 	} else {
# 	  $str=$node;
# 	}
#       } else {
# 	$str=$node;
#       }
#     }
#   else {
#     $seen->{$node}=1;
#     my @sons=$self->{graph}->successors($node);
#     $str="(";
#     foreach my $son (@sons) {
#       $str=$str.$self->eNewick_full_aux($son,$seen,$node).",";
#     }
#     chop($str);
#     if ($self->is_hybrid_node($node)) {
#       my $tag=substr $node,1;
#       if ((defined $parent) &&
# 	  ($self->graph->has_edge_attribute($parent,$node,'type'))) {
# 	$str.=')#'.$self->graph->get_edge_attribute($parent,$node,'type').$tag;
#       } else {
# 	$str.=")$node";
#       }
#     } else {
#       $str.=")$node";
#     }
#   }
#   if ((defined $parent) &&
#       ($self->graph->has_edge_weight($parent,$node))) {
#     $str.=":".$self->graph->get_edge_weight($parent,$node);
#   }
#   return $str;
# }


# displaying data

use overload '""' => \&display;

=head2 display

 Title   : display
 Usage   : my $str=$net->display()
 Function: returns a string containing all the available information on $net
 Returns : string
 Args    : none

=cut

sub display {
  my ($self)=@_;
  my $str="";
  my $graph=$self->{graph};
  my @leaves=$self->leaves();
  my @nodes=@{$self->{nodes}};
  $str.= "Leaves:\t@leaves\n";
  $str.= "Nodes:\t@nodes\n";
  $str.= "Graph:\t$graph\n";
  $str.= "eNewick:\t".$self->eNewick()."\n";
  $str.= "Full eNewick:\t".$self->eNewick_full()."\n";
  $str.= "Mu-data and heights:\n";
  foreach my $node (@nodes) {
    $str.= "v=$node: ";
    if (exists $self->{labels}->{$node}) {
      $str.="\tlabel=".$self->{labels}->{$node}.",";
    } else {
      $str.="\tlabel=(none),";
    }
    $str.= "\th=".$self->{h}->{$node}.", \tmu=".$self->{mudata}->{$node}."\n";
  }
  if (exists $self->{has_temporal_representation}) {
    $str.= "Temporal representation:\n";
    if ($self->{has_temporal_representation}) {
      foreach my $node (@nodes) {
	$str.= "v=$node; ";
	$str.= "\tt=".$self->{temporal_representation}->{$node}."\n";
      }
    } else {
      $str.= "Does not exist.\n";
    }
  }
  if (exists $self->{tripartitions}) {
    $str.= "Tripartitions:\n";
    foreach my $node (@nodes) {
      $str.= "v=$node; ";
      $str.= "\ttheta=".$self->{tripartitions}->{$node}."\n";
    }
  }
  return $str;
}

1;
