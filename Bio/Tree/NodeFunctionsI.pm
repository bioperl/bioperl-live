#
# BioPerl module for Bio::Tree::TreeFunctionsI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org> and Greg Jordan <gjuggler-at-gmail-dot-com>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodeFunctionsI - Decorated Interface implementing basic Node methods

=head1 SYNOPSIS

  # Load a Bio::Tree from a Newick-formatted string.
  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-string => '(a,(b,c));');
  my $tree = $in->next_tree;

  # Get the root node of the tree as a Bio::Tree::Node object
  my $root_node = $tree->root;
  # Human-readable ASCII representation of the whole tree.
  print $root_node->ascii;

  my $node_b = $root_node->find('b');
  $tree->reroot_above($node_b);
  print $tree->find('a')->depth_to_root . "\n";

  # String / file output.
  print $tree->newick."\n";
  $tree->to_file("my_tree.nh");
  
=head1 DESCRIPTION

This interface provides a set of methods providing tree manipulation
and querying functionality, building up only from the five methods
defined by the Bio::Tree::NodeI interface. Any class which implements
the NodeI interface can gain access to the NodeFunctionsI methods by
simply inheriting from NodeFunctionsI as well as NodeI, e.g. 'use base
qw(Bio::Tree::NodeI Bio::Tree::NodeFunctionsI)'



=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich, Aaron Mackey, Justin Reese, Greg Jordan

Email jason-at-bioperl-dot-org
Email amackey-at-virginia.edu
Email jtr4v-at-virginia.edu
Email gjuggler-at-gmail-dot-com

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

Rerooting code was worked on by

  Daniel Barker d.barker-at-reading.ac.uk
  Ramiro Barrantes Ramiro.Barrantes-at-uvm.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tree::NodeFunctionsI;
use strict;


=head2 as_text

 Title   : as_text
 Usage   : my $tree_as_string = $tree->as_text($format)
 Function: Returns the tree as a string representation in the 
           desired format (currently 'newick', 'nhx', or 
           'tabtree')
 Returns : String, the tree as expressed in the desired format
 Args    : String, format type as specified by Bio::TreeIO
 Note    : This method loads the Bio::TreeIO::$format module
           on the fly, and commandeers the _write_tree_Helper
           routine therein to create the tree string. 

=cut

sub as_text {
    my $self = shift;
    my $format = shift;
    my $params_input = shift || {};

    $format = 'nhx' unless (defined $format);

    my $iomod = "Bio::TreeIO::$format";
    $self->_load_module($iomod);

    my $string = '';
    open(my $fh,">",\$string) or die ("Couldn't open $string as file: $!\n");
    my $test = $iomod->new(-format=>$format,-fh=>$fh);

    $test->set_params($params_input);
    $test->write_tree($self);
    close($fh);
    return $string;
}

sub to_newick {
  my $self = shift;

  return $self->as_text('nhx');

}

=head2 translate_ids

 Title   : translate_ids
 Usage   : $node->translate_ids({a => 'Aardvark', b => 'Baobab'});
 Function: Translates the IDs of all nodes beneath the given node based
           on the key-value names of an input hash reference.
 Returns : Nothing (IDs are translated in place)
 Args    : Hashmap, containing key-value pairs where the key is the old
           ID of a node and the value is the desired new ID

=cut
sub translate_ids {
    my $self = shift;
    my $id_map = shift;

    foreach my $node ($self->nodes) {
	my $new_id = $id_map->{$node->id};
	if (defined $new_id) {
	    $node->id($new_id);
	}
    }
}


sub remove_internal_node_labels {
    my $self = shift;

    foreach my $node ($self->nodes) {
	if (!$node->is_leaf) {
	    $node->id('');
	}
    }
    if (!$self->is_leaf) {
	$self->id('');
    }
}

=head2 is_leaf

 Title   : is_leaf
 Usage   : do_something_only_to_leaves() if ($node->is_leaf);
 Function: Return true if this node is a leaf node. This
           implementation counts the number of nodes in the
           node->children() array to determine leaf status.
 Returns : 1 if this node is a leaf, 0 if not
 Args    : None

=cut

sub is_leaf{
    return scalar(shift->children) == 0;
}
sub is_Leaf { shift->is_leaf(@_) }

=head2 depth_to_root

 Title   : depth_to_root
 Usage   : my $depth = $node->depth_to_root();
 Function: Return the number of nodes between this node and the root,
           including this node (i.e. $tree->root->depth_to_root() is 1)
 Returns : Integer
 Args    : None

=cut

sub depth_to_root {
   my ($self) = @_;
   
   my $depth = 1;
   my $node = $self;
   while( defined $node->ancestor ) { 
       $depth += 1;
       $node = $node->ancestor;
   }
   return $depth;
}

=head2 height_to_root

 Title   : height_to_root
 Usage   : my $height = $node->height_to_root();
 Function: Return the total branch length between this node and the
           root node, not including the branch length (if any) above
           the root node
 Returns : Float
 Args    : None

=cut

sub height_to_root {
    my ($self) = @_;

   my $height = 0;
   my $node = $self;
   while( defined $node->ancestor ) { 
       $height += $node->branch_length;
       $node = $node->ancestor;
   }
   return $height;
}
sub distance_to_root { shift->height_to_root(@_) }

=head2 max_distance_to_leaf

 Title   : max_distance_to_leaf
 Usage   : my $max_dist = $node->max_distance_to_leaf();
 Function: Return the maximum total branch length between this node
           and any leaf in the subtree beneath this node
 Returns : Float
 Args    : None

=cut

sub max_distance_to_leaf{
    my ($self) = @_;

    return 0 if( $self->is_leaf );
    
    my $max = 0;
    foreach my $subnode ( $self->children ) { 
	# Collect a 'safe' version of the branch length which returns a 'default' 
	# 1 if the branch length is not defined.
	my $s = $subnode->max_distance_to_leaf + $subnode->branch_length_or_one;
	if( $s > $max ) { $max = $s; }
    }
    return $max;
}
sub max_height_to_leaf { shift->max_distance_to_leaf(@_) }
sub height { shift->max_distance_to_leaf(@_) }

=head2 max_depth_to_leaf

 Title   : max_depth_to_leaf
 Usage   : my $max_depth = $node->max_depth_to_leaf;
 Function: Return the maximum number of nodes between this node and
           any leaf in the subtree beneath this node (including the
           leaf, e.g., $leaf_node->max_depth_to_leaf == 1)
 Returns : Integer
 Args    : None

=cut

sub max_depth_to_leaf{
    my ($self) = @_;

    return 1 if( $self->is_leaf );
    
    my $max = 0;
    foreach my $subnode ( $self->children ) { 
	my $s = $subnode->max_depth_to_leaf;
	if( $s > $max ) { $max = $s; }
    }
    return $max;
}
sub depth { shift->max_depth_to_leaf(@_) }

=head2 total_branch_length

 Title   : total_branch_length
 Usage   : my $total_length = $node->total_branch_length;
 Function: Return the total summed branch length of the subtree
           beneath the current node, NOT including the length of the
           branch directly above the current node.
 Returns : Float
 Args    : None

=cut

sub total_branch_length {
    my $self = shift;
    my $sum = 0;
    foreach my $node ( $self->nodes ) {
      $sum += $node->branch_length || 0;
    }
    return $sum;
}
sub subtree_length { shift->total_branch_length(@_) }
sub subtree_size { shift->total_branch_length(@_) }

sub children_branch_length {
  my $self = shift;

  my $bl = 0;
  my @children = $self->children;
  foreach my $child (@children) {
    $bl += $child->branch_length;
  }

  return $bl;
}

=head2 distance

 Title   : distance
 Usage   : my $dist = $node_a->distance($node_b);
 Function: Return the branch length distance between this node and another node in the tree.
 Returns : Float, the distance between the two nodes
 Args    : Bio::Tree::NodeI, the other node

=cut

sub distance {
    my $self = shift;
    my $other_node = shift;

    my $lca = $self->lca($other_node);

    my $lca_distance = $lca->distance_to_root;
    my $my_distance = $self->distance_to_root;
    my $other_distance = $other_node->distance_to_root;

    return ($my_distance - $lca_distance) + ($other_distance - $lca_distance);
}

=head2 lca

 Title   : lca
 Usage   : my $lca = $node->lca($other_node);
 Function: Return the lowest common ancestor (aka most recent common
           ancestor) between this node and another node in the tree.
 Returns : Bio::Tree::NodeI, the lowest common ancestor
 Args    : Bio::Tree::NodeI, the other node

=cut

sub lca {
    my $self = shift;
    my @other_nodes = @_;

    my $other_node = shift @other_nodes;
    foreach my $other_other (@other_nodes) {
      $other_node = $other_node->lca($other_other);
    }

    my @lineage = $self->lineage;
    my @other_lineage = $other_node->lineage;

    my $seen_ancestors;
    my $shared_ancestors;

    # Include $self and $other_node in the LCA calculation, as one
    # node might actually be an ancestor of the other

    foreach my $node ($self,$other_node,@lineage,@other_lineage) {
	if (!defined $seen_ancestors->{$node}) {
	    $seen_ancestors->{$node} = $node;
	} else {
	    # We've already seen this ancestor once, which means it's a shared ancestor.
	    # Calculate and store its depth to the root node.
	    $shared_ancestors->{$node} = $node;
	}
    }

    # Every node within shared_ancestor_depths is shared, now sort by
    # depth to find the deepest, aka lowest, common ancestor.
    my $max_depth = 0;
    my $max_node;

    foreach my $key (keys %$shared_ancestors) {
	my $node = $shared_ancestors->{$key};
	my $depth = $node->depth_to_root;
	if ($depth > $max_depth) {
	    $max_depth = $depth;
	    $max_node = $node;
	}
    }

    return $max_node;
}

=head2 leaves

 Title   : leaves
 Usage   : my @leaves = $node->leaves;
 Function: Return a list of all the leaf nodes contained in the
           subtree beneath the current node. Returns a list containing
           only the current node if called on a leaf node.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub leaves {
    my $self = shift;

    my @all_nodes = $self->nodes();
    return grep {$_->is_leaf} @all_nodes;
}
sub get_all_leaves { shift->leaves(@_) }

=head2 nodes

 Title   : nodes
 Usage   : my @nodes = $node->nodes;
 Function: Return a list of all nodes contained in the subtree beneath
           the current node. Returns a list containing only the
           current node if called on a leaf node.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub nodes {
   my ($self) = @_;

   my @nodes = ($self);
   foreach my $node ( $self->children() ) {
       push @nodes, ($node->nodes());
   }
   return @nodes;
}
sub get_all_Descendents { shift->nodes(@_) }
sub get_Descendents { 
    my $self = shift;
    $self->deprecated(-message => 'get_Descendents() is ambiguous and deprecated... use nodes() instead');
    $self->nodes(@_);
}
sub get_all_nodes { shift->nodes(@_) }

=head2 _ordered_nodes

 Title   : _ordered_nodes
 Usage   : my @nodes = $node->_ordered_nodes('depth')
 Function: Helper method for nodes_breadth_first() and
           nodes_depth_first(). Returns a list of all nodes in the
           subtree beneath the current node, ordered by either
           breadth-first or depth-first traversal order.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : 'depth' or 'breadth' to indicate desired traversal order (defaults to 'depth')

=cut

sub _ordered_nodes {
    my $self = shift;
    my $order = shift;
   
    $order ||= 'depth';

    if ($order =~ m/^b|(breadth)$/oi) {
       my @children = ($self);
       for (@children) {
        push @children, $_->children;
       }
       return @children;
   }

   if ($order =~ m/^d|(depth)$/oi) {
       # this is depth-first search I believe
       my @children = ($self->nodes);
       return @children;
   }
}

=head2 nodes_breadth_first

 Title   : nodes_breadth_first
 Usage   : my @nodes = $node->nodes_breadth_first();
 Function: Return a list of all nodes in the subtree beneath the
           current node, in breadth-first traversal order.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub nodes_breadth_first {
    my $self = shift;
    return $self->_ordered_nodes('b');
}

=head2 nodes_depth_first

 Title   : nodes_depth_first
 Usage   : my @nodes = $node->nodes_depth_first();
 Function: Return a list of all nodes in the subtree beneath the
           current node, in depth-first traversal order.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub nodes_depth_first {
    my $self = shift;
    return $self->_ordered_nodes('d');
}

=head2 lineage

 Title   : lineage
 Usage   : my @lineage_from_root = $node->lineage();
 Function: Return a list of all nodes in the lineage between the
           current node and the root node, starting with the root node
           and NOT including the current node.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub lineage {
    my $self = shift;

    my @lineage = ();
    my $node = $self->parent;
    while ($node) {
        unshift(@lineage, $node);
        $node = $node->parent;
    }
    return @lineage;
}

=head2 change_child_to_parent

 Title   : change_child_to_parent
 Usage   : $parent_node->change_child_to_parent($child_node);
 Function: Turns the indicated child of this node into the parent of
           this node, bubbling this edge direction-flipping up the
           tree all the way to the root node. This method implements
           the main transformation involved in the Bio::Tree::Tree
           reroot() method.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the child node which will become the new
           parent of this node

=cut

sub change_child_to_parent {
    my $self = shift;
    my $child_to_become_parent = shift;

    if ($child_to_become_parent->parent != $self) {
	$self->warn("Called 'change_child_to_parent' using a non-child node -- check your input! Doing nothing...");
	return;
    }

    if ($self->parent) {
	# Climb up the tree, flipping the direction of all nodes up to the root.
	$self->parent->change_child_to_parent($self);
    }

    # Store the branch length between us and our child.
    my $bl = $child_to_become_parent->branch_length;

    # Remove $child from me, add me to $child
    $self->remove_child($child_to_become_parent);
    $child_to_become_parent->add_child($self);
    # Set my branch length appropriately
    $self->branch_length($bl);
}

=head2 child_count

 Title   : child_count
 Usage   : my $num_direct_children = $node->child_count;
 Function: Return the number of children directly beneath the current node
 Returns : Integer, the number of direct children.
 Args    : None

=cut

sub child_count {
    my $self = shift;
    return scalar($self->children);
}

=head2 node_count

 Title   : node_count
 Usage   : my $number_of_nodes_beneath = $node->node_count;
 Function: Return the total number of nodes in the subtree beneath the
           current node, not including the curent node itself.
 Returns : Integer, the number of nodes in the subtree beneath the
           current node
 Args    : None

=cut

sub node_count {
   my ($self) = @_;
   return scalar($self->nodes);
}
sub descendent_count { shift->node_count(@_) }
sub number_nodes { shift->node_count(@_) }

=head2 leaf_count

 Title   : leaf_count
 Usage   : my $number_of_leaves_beneath = $node->leaf_count;
 Function: Return the total number of leaf nodes in the subtree
           beneath the current node, not including the curent node
           itself.
 Returns : Integer, the number of leaf nodes in the subtree beneath
           the current node
 Args    : None

=cut

sub leaf_count {
   my ($self) = @_;
   return scalar($self->leaves);
}

=head2 split_branch_with_node

 Title   : split_branch_with_node
 Usage   : $node->split_branch_with_node($new_node,0.5);
 Function: Split the branch above the current node at a certain
           position by splicing the passed node into the tree, Ex:
           given (A---->B) and calling
           B->split_branch_with_node(C, 0.5), we end up with
           (A-->C-->B).
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the node to insert into the branch above the
           current node
           Float, the fractional distance along the branch above the
           current node at which to place the new node

=cut

sub split_branch_with_node {
    # Splits the branch above this node with the given node.
    my $self = shift;
    my $node_to_insert = shift;
    my $fraction_above_node = shift;

    $fraction_above_node = 0.5 unless (defined $fraction_above_node);

    if (!defined $self->parent) {
	$self->warn("Cannot split branch above the root node -- Not inserting node!");
	return undef;
    }

    my $parent = $self->parent;
    my $bl = $self->branch_length;

    $parent->add_child($node_to_insert);
    $parent->remove_child($self);
    $node_to_insert->add_child($self);

    $node_to_insert->branch_length($bl * (1-$fraction_above_node));
    $self->branch_length($bl * $fraction_above_node);
    
    return $node_to_insert;
}

=head2 remove_children

 Title   : remove_children
 Usage   : $node->remove_children;
 Function: Remove all of the children directly beneath this node from
           the tree
 Returns : Nothing
 Args    : None

=cut

sub remove_children{
   my ($self) = @_;

   foreach my $node ($self->children) {
       $self->remove_child($node);
   }   
}
sub remove_Descendents { shift->remove_children(@_) }

=head2 remove

 Title   : remove
 Usage   : $node->remove();
 Function: Removes this node and its contained subtree from the tree
           by calling remove_child on this node's parent
 Returns : Nothing
 Args    : None

=cut

sub remove {
    my $self = shift;
    if ($self->parent) {
	$self->parent->remove_child($self);
    }
}
sub remove_from_parent { shift->remove(@_) }

=head2 splice

 Title   : splice
 Usage   : $node->splice();
 Function: Splices the current node out of the tree by removing this
           node and re-attaching its former children to its former
           parent, preserving the lineage branch length for each node.
 Returns : Nothing
 Args    : None

=cut

sub splice {
    my $self = shift;

    if ($self->parent) {
	my $parent = $self->parent;
	my $bl = $self->branch_length;
	foreach my $child ($self->children) {
	    # Remove the child from me, add it to my parent.
	    $self->remove_child($child);
	    $parent->add_child($child);

	    # Add my branch length to my child -- hence this method is called "splice"
	    my $child_bl = $child->branch_length;
	    if (defined $bl && defined $child_bl) {
		$child->branch_length($child_bl + $bl);
	    } elsif (defined $child_bl) {
		$child->branch_length($child_bl);
	    } elsif (defined $bl) {
		$child->branch_length($bl);
	    }
	}
	$parent->remove_child($self);
    } else {
      warn("Called splice() on a node with no parent -- sure it's not the root node?");
    }
}

=head2 force_binary

 Title   : force_binary
 Usage   : $node->force_binary; print $node->ascii;
 Function: Break all multifurcations in the subtree below this node by
           separating out each multifurcation into a series of
           bifurcations
 Returns : Nothing
 Args    : None

=cut

sub force_binary {
    my $self = shift;

    foreach my $child ($self->children) {
	# Recurse through tree.
	$child->force_binary;
    }

    my @children = $self->children;
    while (scalar(@children) > 2) {
	my $first_child = shift @children;
	
	# Make a new node to hold a binary split.
	my $new_node = new $self;
	$first_child->parent->add_child($new_node);
	foreach my $cur_child (@children) {
	    # Remove the "extra" children from the direct parent
	    $cur_child->remove;
	    # and add them to $new_node
	    $new_node->add_child($cur_child);
	}
    }
}

=head2 is_binary

 Title   : is_binary
 Usage   : print "Node is ok!\n" if ($node->is_binary);
 Function: Return true if the current node is binary (without
           considering sub-trees). Call is_subtree_binary() to test
           for a completely binary subtree.
 Returns : True or false
 Args    : None

=cut

sub is_binary {
    my $self = shift;
    return ($self->is_leaf || scalar($self->children) == 2);
}

=head2 is_subtree_binary

 Title   : is_subtree_binary
 Usage   : print "Subtree OK!\n" if ($node->is_subtree_binary);
 Function: Return true if the entire subtree beneath the current node
           is binary.
 Returns : True or false
 Args    : None

=cut

sub is_subtree_binary {
    my $self = shift;

    foreach my $node ($self->nodes) {
	return 0 unless ($node->is_binary);
    }
    return 1;
}

=head2 contract_linear_paths

 Title   : contract_linear_paths
 Usage   : $node->contract_linear_paths
 Function: Splice out linear nodes (nodes with exactly 1 parent and 1
           child) from the tree. The root node is never considered
           linear, so it won't be affected by this method.
 Returns : Nothing
 Args    : None

=cut

sub contract_linear_paths {
    my $self = shift;
    my $preserve_root = shift;

    if (!defined $preserve_root) {
      $preserve_root = 0;
    }

    foreach my $child ($self->children) {
	$child->contract_linear_paths;
    }

    if ($self->child_count == 1) {
      if ($preserve_root && !$self->parent) {
        # Do nothing -- we're a root node, but want to be kept!
      } elsif (!$self->parent) {
        # We're a root node with an only child -- splice out the child...
        my @children = $self->children;
        my $only_child = $children[0];
        $self->branch_length($only_child->branch_length);
        $only_child->branch_length(0);
        $only_child->splice;
      } else {
        $self->splice;
      }
    }
}
sub remove_elbow_nodes { shift->contract_linear_paths(@_) }

=head2 branch_length_or_one

 Title   : branch_length_or_one
 Usage   : my $bl = $node->branch_length_or_one;
 Function: Return the node's branch length if defined, or 1 if
           undefined. Useful for situations where trees without
           explicit branch lengths should be assumed to have all
           1-length branches.
 Returns : Nothing
 Args    : None

=cut


sub branch_length_or_one {
    my $self = shift;
    
    my $bl = $self->branch_length;
    return 1 unless (defined $bl);
    return $bl;
}

=head2 ascii
 Title   : ascii
 Usage   : print STDERR $tree->ascii."\n";
 Function: Return a representation of this tree as ASCII text. Shamelessly
           copied / ported from PyCogent code. http://pycogent.sourceforge.net/
 Returns : A multi-line String, suitable for printing to the console
 Args    : show_internal - 0 to hide internal nodes, 1 to show them (default 1)
           compact - 1 to use a 'compact' mode with exactly 1 line per node. (default 0)
           ignore_branchlengths - 0 to ignore branch lengths, 1 to use them. (default 0)
=cut

sub ascii {
    my $self = shift;
    my $show_internal = shift;
    my $compact = shift;
    my $ignore_branchlengths = shift;

    my $max_bl = $self->max_distance_to_leaf;
    if ($max_bl == 0) {
	$max_bl = $self->max_depth_to_leaf;
    }
    my $width = 80;

    my $char_per_bl = int($width / $max_bl);

    my ($lines_arrayref,$mid) = $self->_ascii_art($self,'-',$char_per_bl,$show_internal,$compact,$ignore_branchlengths);
    my @lines = @{$lines_arrayref};
    return join("\n",@lines)."\n";
}

# Private helper method for $tree->ascii_art
sub _ascii_art {
    my $self = shift;
    my $node = shift;
    my $char1 = shift;
    my $char_per_bl = shift;
    my $show_internal = shift;
    my $compact = shift;
    my $ignore_branchlengths = shift;

    $char1 = '-' unless (defined $char1);
    $show_internal = 1 unless (defined $show_internal);
    $compact = 0 unless (defined $compact);
    $ignore_branchlengths = 1 unless (defined $ignore_branchlengths);

    my $len = 10;
    if (!$ignore_branchlengths) {
	my $bl = $node->branch_length;
	$bl = 1 unless (defined $bl && $bl =~ /^\-?\d+(\.\d+)?$/);
	$len = int($bl * $char_per_bl);
	$len = 1 if ($len < 1);
    }
    my $pad = ' ' x $len;
    my $pa = ' ' x ($len-1);
    my $name_str = $node;
    if ($node->isa("Bio::Tree::LabeledNodeI")) {
	$name_str = $node->id;
    }
    $name_str = '' unless (defined $name_str);

    if (!$node->is_leaf) {
	my @mids;
	my @results;
	my @children = $node->children;
	for (my $i=0; $i < scalar(@children); $i++) { 
	    my $child = $children[$i];
	    my $char2;
	    if ($i == 0) {
		# First child.
		$char2 = '/';
	    } elsif ($i == scalar(@children) -1) {
		# Last child.
		$char2 = '\\';
	    } else {
		# Middle child.
		$char2 = '-';
	    }
	    my ($clines_arrayref, $mid) = $self->_ascii_art($child,$char2,$char_per_bl,$show_internal,$compact,$ignore_branchlengths);
	    push @mids,($mid+scalar(@results));
	    my @clines = @$clines_arrayref;
	    foreach my $line (@clines) {
		push @results, $line;
	    }
	    if (!$compact) {
		push @results,'';
	    }
	}
	if (!$compact) {
	    pop @results;
	}
	my $lo = $mids[0];
	my $hi = $mids[scalar(@mids)-1];
	my $end = scalar(@results);

	my @prefixes = ();
	push @prefixes, ($pad) x ($lo+1);
	push @prefixes, ($pa.'|') x ($hi-$lo-1);
	push @prefixes, ($pad) x ($end-$hi);
	
	my $mid = int(($lo + $hi) / 2);
	$prefixes[$mid] = $char1 . '-'x($len-2) . substr($prefixes[$mid],length($prefixes[$mid])-1,1);
	my @new_results;
	for (my $i=0; $i < scalar(@prefixes); $i++) {
	    my $p = $prefixes[$i] || '';
	    my $r = $results[$i] || '';
	    push @new_results, ($p . $r);
	}
	@results = @new_results;
	if ($show_internal) {
	    my $stem = $results[$mid];
	    my $str = substr($stem,0,1);
	    $str .= $name_str;
	    my $substr_index = length($name_str)+1;
	    $substr_index = length($stem) if ($substr_index > length($stem));
	    my $stem_substr = substr($stem,$substr_index);
	    $stem_substr = '' unless (defined $stem_substr);
	    $str .= $stem_substr;
	    $results[$mid] = $str;
	}
	return (\@results,$mid);
    } else {
	my @results = ($char1 . ('-'x$len) . $name_str);
	if ($ignore_branchlengths) {
	    @results = ($char1 . '-' . $name_str);
	}
	return (\@results,0);
    }
}

=head2 find_by_id

 Title   : find_by_id
 Usage   : my $node_a = $root_node->find_by_id('a');
 Function: Search through the tree for a given ID, returning all matching nodes.
 Returns : (array context) An array of all nodes with a matching ID
           (scalar context) the first node with the matching ID
 Args    : String, the id for which to search through the tree

=cut

sub find_by_id {
    my $self = shift;
    my $value = shift;

    my @nodes;
    foreach my $node ($self,$self->nodes) {
	push @nodes, $node if ($node->isa("Bio::Tree::LabeledNodeI") && $node->id eq $value);
    }
    
    if ( wantarray) { 
	return @nodes;
    } else { 
	if( @nodes > 1 ) { 
	    $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
	return shift @nodes;
    }
}
sub find { shift->find_by_id(@_) }

=head2 find_by_tag_value

 Title   : find_by_tag_value
 Usage   : my $node_a = $root_node->find_by_tag_value('a');
 Function: Search through the tree for a given tag value, returning all matching nodes.
 Returns : (array context) An array of all nodes with a matching tag value
           (scalar context) the first node with the matching tag value
 Args    : String, the tag to search within
           String, the value to search for

=cut

sub find_by_tag_value {
    my $self = shift;
    my $tag = shift;
    my $value = shift;

    my @nodes;
    foreach my $node ($self->nodes) {
	next unless ($node->isa("Bio::Tree::TagValueHolder"));
	my @values = $node->get_tag_values($tag);
	if (grep {$_ eq $value} @values) {
	    push @nodes, $node ;
	}
    }

    if ( wantarray) { 
	return @nodes;
    } else { 
	if( @nodes > 1 ) { 
	    $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
	return shift @nodes;
    }
}

=head2 find_by_tag_regex

 Title   : find_by_tag_regex
 Usage   : my $node_a = $root_node->find_by_tag_value('species', '.*coli.*');
 Function: Search through the tree for a given tag value using a regular expression,
           returning all nodes whose tag value matches the regex.
 Returns : (array context) An array of all nodes with a matching tag value
           (scalar context) the first node with the matching tag value
 Args    : String, the tag to search within
           String, the regular expression to search for

=cut

sub find_by_tag_regex {
    my $self = shift;
    my $tag = shift;
    my $value_regex = shift;

    my @nodes;
    foreach my $node ($self->nodes) {
	next unless ($node->isa("Bio::Tree::TagValueHolder"));
	my @values = $node->get_tag_values($tag);
	if (grep {$_ =~ $value_regex} @values) {
	    push @nodes, $node ;
	}
    }

    if ( wantarray) { 
	return @nodes;
    } else { 
	if( @nodes > 1 ) { 
	    $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
	return shift @nodes;
    }
}


sub enclosed_leaves_string {
  my $self = shift;
  my $sep = shift;

  $sep = '|' unless (defined $sep);

  if ($self->is_leaf) {
    return $self->id;
  }

  my @leaves_beneath = map {$_->id} $self->leaves;
  my $leaf_string = join($sep, sort {$a cmp $b} @leaves_beneath);
  return $leaf_string;
}


# Extracts the minimum spanning subtree defined by the given nodes. Nodes can either be 
# leaf or internal nodes, and by default all 'elbow' nodes left after the subtree extraction
# are spliced out of the tree. To keep the internal nodes, call the method
# subtree_with_internals
# Note: this CLONES the tree and returns a slice of the clone, so the original tree structure
# is not modified!
sub slice {
    my $self = shift;
    return $self->_slice(1, @_);
}

sub slice_by_ids {
    my $self = shift;
    my @ids = @_;

    my @nodes;
    foreach my $id (@ids) {
	my $node = $self->find($id);
        if ($node) {
          push @nodes, $node;
        } else {
          warn("Node $id not found during tree slice operation!");
        }
    }
    return $self->slice(@nodes);
}

sub slice_with_internals {
    my $self = shift;
    return $self->_slice(0, @_);
}

sub _get_equivalent_nodes_in_cloned_tree {
    my $self = shift;
    my $cloned_tree = shift;
    my @node_list = @_;

    my @node_indices;
    my @nodes = $self->nodes;
    for (my $i=0; $i < scalar(@nodes); $i++) {
	my $cur_node = $nodes[$i];
	if (grep {$cur_node == $_} @node_list) {
	    push @node_indices, $i;
	}
    }
    
    my @cloned_list;
    my @cloned_nodes = $cloned_tree->nodes;
    foreach my $index (@node_indices) {
	push @cloned_list, $cloned_nodes[$index];
    }
    return @cloned_list;
}

sub _slice {
    my $self = shift;
    my $remove_internals = shift;
    my @nodes_to_keep = @_;

    # Create a clone of the current tree, and match up the original nodes_to_keep
    # with their equivalent in the new tree.
    my $clone = $self->clone;
    @nodes_to_keep = $self->_get_equivalent_nodes_in_cloned_tree($clone, @nodes_to_keep);

    # Keep a hash where the keys are refs to the nodes we want to keep.
    my $keepers_hash;
    foreach my $node (@nodes_to_keep) {
	$keepers_hash->{$node} = 1;
	# Get the array of all nodes from the root to this node.
	my @lineage = $node->lineage;
	foreach my $lineage_node (@lineage) {
	    $keepers_hash->{$lineage_node} = 1;
	}
    }
    
    # Create an array of nodes to remove. Don't remove a node if it's in the
    # keepers hash.
    my @remove_me;
    foreach my $node ($clone->nodes) {
	if (!defined $keepers_hash->{$node}) {
	    push @remove_me, $node;
	}
    }

    # Go through each node in the remove_me list and remove it from the tree.
    foreach my $remove_node (@remove_me) {
	$remove_node->remove;
    }

    if ($remove_internals) {
	$clone->contract_linear_paths;
        $clone->branch_length(0);
    }
    return $clone;
}

1;
