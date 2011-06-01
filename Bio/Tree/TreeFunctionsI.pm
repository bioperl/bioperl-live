#
# BioPerl module for Bio::Tree::TreeFunctionsI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeFunctionsI - Decorated Interface implementing basic Tree-wide methods

=head1 SYNOPSIS


=head1 DESCRIPTION

This interface provides a set of implementated tree manipulation
functions which depend only on the methods defined in the the TreeI and NodeI
interfaces.

Many methods which perform tree operations simply forward the method
call to the root node, meaning that most tree operations are actually
implemented in the NodeFunctionsI module, even if they are accessible
to a Bio::Tree::Tree object. Only methods which require manipulation
of the link between the Bio::Tree::TreeI object and the root
Bio::Tree::NodeI, such as reroot, should be implemented here.

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

=head1 AUTHOR - Jason Stajich, Aaron Mackey, Justin Reese

Email jason-at-bioperl-dot-org
Email amackey-at-virginia.edu
Email jtr4v-at-virginia.edu

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

package Bio::Tree::TreeFunctionsI;
use strict;

sub slice { 
  return new Bio::Tree::Tree(-root => shift->root->slice(@_));
}
sub slice_by_ids { 
  return new Bio::Tree::Tree(-root => shift->root->slice_by_ids(@_));
}

sub print_tree { print shift->ascii(@_) }
sub ascii { shift->root->ascii(@_) }
sub as_text { shift->root->as_text(@_) }
sub to_string { shift->root->as_text(@_) }
sub to_newick { shift->root->to_newick(@_) }
sub newick { shift->root->to_newick(@_) }

sub translate_ids { shift->root->translate_ids(@_) }

sub remove_internal_node_labels{ shift->root->remove_internal_node_labels(@_) }

sub contract_linear_paths { shift->root->contract_linear_paths(@_) }

sub get_lineage_nodes {
    my $self = shift;
    my $node = shift;
    return $node->lineage;
}

sub total_branch_length { shift->root->total_branch_length(@_) }
sub subtree_length { shift->root->subtree_length(@_) }

sub max_distance_to_leaf { shift->root->max_distance_to_leaf(@_) }
sub max_depth_to_leaf { shift->root->max_depth_to_leaf(@_) }

sub number_nodes { shift->root->node_count(@_) } # This alias sucks...

=head2 nodes

 Title   : nodes
 Usage   : my @all_nodes = $tree->nodes
 Function: Returns all nodes within the tree, including the root node
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub nodes {
  my $self = shift;
  return () unless ($self->root);
  return $self->root->nodes;
}
sub get_nodes { shift->nodes }
sub leaf_nodes { shift->root->leaves }
sub leaves { shift->root->leaves }
sub get_leaves { shift->root->leaves }
sub get_leaf_nodes { shift->root->leaves }

sub nodes_breadth_first { shift->root->nodes_breadth_first(@_) }
sub nodes_depth_first { shift->root->nodes_depth_first(@_) }

=head2 remove_node

 Title   : remove_node
 Usage   : $tree->remove_node($node); $tree->remove_node('a');
 Function: If given a Bio::Tree::NodeI, removes that node from the
           tree by calling $node->remove. If given a string, searches
           for the first node with that ID and removes it from the
           tree.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : None

=cut

sub remove_node {
    my $self = shift;
    my $node = shift;
    if (ref $node && $node->isa("Bio::Tree::NodeI")) {
	$node->remove;
    } else {
	my $found_node = $self->find_node($node);
	if ($found_node) {
	    $found_node->remove;
	} else {
	    $self->warn("Could not find node to remove: [$node]\n");
	}
    }
}
sub remove_Node { shift->remove_node(@_) }


sub find { shift->root->find_by_id(@_) }
sub find_node { shift->root->find_by_id(@_) }
sub find_by_id { shift->root->find_by_id(@_) }

sub force_binary { shift->root->force_binary(@_) }
sub is_binary { shift->root->is_subtree_binary(@_) }


=head2 get_lca

 Title   : get_lca
 Usage   : my $lca_node = get_lca(@nodes);
 Function: given two or more nodes, returns the lowest common ancestor (aka most
           recent common ancestor)
 Returns : node object or undef if there is no common ancestor
 Args    : -nodes => arrayref of nodes to test, OR
           just a list of nodes

=cut

sub get_lca {
    my $self = shift;
    my @nodes = @_;
    
    @nodes >= 2 or $self->throw("At least 2 nodes are required");

    # Go through all the nodes in the argument array, keeping track of
    # the current LCA as we go.
    my $lca_node = shift @nodes;
    while (scalar(@nodes) > 0) {
	my $other_node = shift @nodes;
	$lca_node = $other_node->lca($lca_node);
    }

    return $lca_node;
}
sub lca { shift->get_lca(@_) }


=head2 merge_lineage

 Title   : merge_lineage
 Usage   : merge_lineage($node)
 Function: Merge a lineage of nodes with this tree.
 Returns : n/a
 Args    : Bio::Tree::TreeI with only one leaf, OR
           Bio::Tree::NodeI which has an ancestor

 For example, if we are the tree $tree:

 +---B
 |
 A
 |
 +---C

 and we want to merge the lineage $other_tree:

 A---C---D

 After calling $tree->merge_lineage($other_tree), $tree looks like:

 +---B
 |
 A
 |
 +---C---D

=cut

sub merge_lineage {
    my ($self, $thing) = @_;
    $self->throw("Must supply an object reference") unless ref($thing);

    my ($lineage_tree, $lineage_leaf);
    if ($thing->isa('Bio::Tree::TreeI')) {
        my @leaves = $thing->get_leaf_nodes;
        $self->throw("The supplied Tree can only have one leaf") unless @leaves == 1;
        $lineage_tree = $thing;
        $lineage_leaf = shift(@leaves);
    }
    elsif ($thing->isa('Bio::Tree::NodeI')) {
        $self->throw("The supplied Node must have an ancestor") unless $thing->ancestor;
        $lineage_tree = $self->new(-node => $thing);
        $lineage_leaf = $thing;
    }

    # see if any node in the supplied lineage is in our tree - that will be
    # our lca and we can merge at the node below
    my @lineage = ($lineage_leaf, reverse($self->get_lineage_nodes($lineage_leaf)));
    my $merged = 0;
    for my $i (0..$#lineage) {
        my $lca = $self->find_node(-internal_id => $lineage[$i]->internal_id) || next;

        if ($i == 0) {
            # the supplied thing to merge is already in the tree, nothing to do
            return;
        }
        # $i is the lca, so the previous node is new to the tree and should
        # be merged on
        $lca->add_Descendent($lineage[$i-1]);
        $merged = 1;
        last;
    }
    $merged || ($self->warn("Couldn't merge the lineage of ".$lineage_leaf->id." with the rest of the tree!\n") && return);
}


# alias
sub _clone { shift->clone(@_) }

# safe node clone that doesn't seg fault, but deliberately loses ancestors and
# descendents
sub _clone_node {
    my ($self, $node) = @_;
    my $clone = $node->new;

    while (my ($key, $val) = each %{$node}) {
        if ($key eq '_desc' || $key eq '_ancestor') {
            next;
        }
        ${$clone}{$key} = $val;
    }

    return $clone;
}

=head2 distance

 Title   : distance
 Usage   : my $dist = $tree->distance($node_a,$node_b);
 Function: Return the branch length distance between two given nodes in the tree
 Returns : Float, the distance between two nodes
 Args    : Bio::Tree::NodeI, the first node
           Bio::Tree::NodeI, the second node

=cut

sub distance {
    my $self = shift;
    my $node_a = shift;
    my $node_b = shift;
    
    return $node_a->distance($node_b);
}

=head2 is_monophyletic

 Title   : is_monophyletic
 Usage   : die unless ($tree->is_monophyletic([$node_a,$node_b],$outgroup));
 Function: Return true if the given nodes are monophyletic relative to
           the given outgroup. A monophyletic group is defined as a group which
           contains its LCA and all of the LCA's descendants.
 Returns : 1 if the nodes are monophyletic, 0 if monophyly is violated
 Args    : Arrayref of Bio::Tree::NodeI, the set of nodes to test for monophyly
           Bio::Tree::NodeI, the outgroup node

=cut

sub is_monophyletic{
   my $self = shift;
   my $node_arrayref = shift;
   my $outgroup = shift;

   if( ! defined $node_arrayref || ! defined $outgroup ) {
       $self->warn("Must supply nodes and outgroup parameters to the method
is_monophyletic");
       return;
   }
   if( ref($node_arrayref) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid array reference for -nodes");
   }

   my @nodes = @{$node_arrayref};

   my $clade_root = $self->get_lca(@nodes);
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return;
   }

   # Monophyly is violated if an ancestor of the outgroup node is the clade root.
   my $og_ancestor = $outgroup->ancestor;
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor == $clade_root ) {
           # monophyly is violated
           return 0;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }

   # Monophyly is violated if the set of nodes does not contain ALL
   # leaves underneath the clade root.
   my $given_leaves;
   map {$given_leaves->{$_} = $_} @nodes;
   foreach my $leaf ($clade_root->leaves) {
       return 0 unless ($given_leaves->{$leaf});
   }

   return 1;
}

=head2 is_paraphyletic

 Title   : is_paraphyletic
 Usage   : die unless ($tree->is_paraphyletic([$node_a,$node_b],$outgroup));
 Function: Return true if the given set of nodes is paraphyletic. A
           paraphyletic group is an 'incomplete' monophyletic group,
           i.e. a group with one or more descendent species removed
 Returns : 1 if the nodes are paraphyletic, 0 if not
 Args    : Arrayref of Bio::Tree::NodeI, the set of nodes to test for paraphyly
           Bio::Tree::NodeI, the outgroup node

=cut

sub is_paraphyletic{
   my $self = shift;
   my $node_arrayref = shift;
   my $outgroup = shift;

   if( ! defined $node_arrayref || ! defined $outgroup ) {
       $self->warn("Must supply nodes and outgroup parameters to the method is_paraphyletic");
       return;
   }
   if( ref($node_arrayref) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid array reference for -nodes");
       return;
   }

   my @nodes = @{$node_arrayref};

   if ($self->is_polyphyletic(\@nodes,$outgroup)) {
       return 0;
   }

   my $clade_root = $self->get_lca(@nodes);
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return 0;
   }

   # Paraphyletic if the set of nodes does not contain ALL
   # leaves underneath the clade root.
   my $given_leaves;
   map {$given_leaves->{$_} = $_} @nodes;
   foreach my $leaf ($outgroup->parent->leaves) {
       return 1 if (!defined $given_leaves->{$leaf});
   }
   
   return 0;
}

=head2 is_polyphyletic

 Title   : is_polyphyletic
 Usage   : die if ($tree->is_polyphyletic([$node_a,$node_b],$outgroup));
 Function: Return true if the given set of nodes is polyphyletic. A
           polyphyletic group is one with multiple ancestries with
           respect to the given outgroup (e.g., the LCA of the given
           species is ancestral to the outgroup node)
 Returns : 1 if the nodes are polyphyletic, 0 if not
 Args    : Arrayref of Bio::Tree::NodeI, the nodes to test for polyphyly
           Bio::Tree::NodeI, the outgroup node

=cut

sub is_polyphyletic {
   my $self = shift;
   my $node_arrayref = shift;
   my $outgroup = shift;
   my @nodes = @{$node_arrayref};
   my $clade_root = $self->get_lca(@nodes);

   # Polyphyletic if an ancestor of the outgroup node is the clade root.
   my $og_ancestor = $outgroup->parent;
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor == $clade_root ) {
	   return 1;
       }
       $og_ancestor = $og_ancestor->parent;
   }
   return 0;
}

=head2 reroot

 Title   : reroot
 Usage   : $tree->reroot($new_root_node);
 Function: Re-roots the tree so the passed node becomes the new root
           node. No restrictions are placed on the allowable types of
           node for the new root, so you may reroot the tree on a leaf
           node (although that wouldn't make much biological
           sense). Most users probably want to use reroot_above to
           create a new root node along the branch above a given node,
           causing the given node to be one of the two children of the
           new root.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the node to become root

=cut

sub reroot {
    my ($self,$new_root) = @_;
    unless (defined $new_root && $new_root->isa("Bio::Tree::NodeI")) {
        $self->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
        return 0;
    }

    my $old_root = $self->get_root_node;
    if( $new_root == $old_root || !$new_root->ancestor ) {
	$self->warn("Node requested for reroot is already the root node!");
	return 0;
    }
    
    # Note -- this implementation neither adds nor removes any nodes!
    # Most practical re-rootings will require adding 

    $new_root->parent->change_child_to_parent($new_root);
    $new_root->branch_length(undef);
    $self->set_root_node($new_root);

    return 1;
}

=head2 reroot_above

 Title   : reroot_above
 Usage   : $tree->reroot_above($outgroup,0.5);
 Function: Re-roots the tree by creating a new node along the branch
           directly above the given node and re-rooting the tree along this
           (newly created) node.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the node above which to reroot the tree.
           Float, the fractional branch length towards the parent to
           insert the new node (0 inserts at the passed node, 1
           inserts at its parent)

=cut


sub reroot_above {
    my $self = shift;
    my $node = shift;
    my $fraction_above_node = shift;

    $fraction_above_node = 0.5 unless (defined $fraction_above_node);

    unless (defined $node && $node->isa("Bio::Tree::NodeI")) {
        $self->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
        return 0;
    }
    my $old_root = $self->get_root_node;
    if( $node == $old_root || !$node->ancestor ) {
	$self->warn("Node requested for reroot is already the root node!");
	return 0;
    }

    # Reroot the tree at a new node located halfway on branch above $node
    my $new_root = new $node;
    $node->split_branch_with_node($new_root,$fraction_above_node);

    return $self->reroot($new_root);
}

=head2 move_id_to_bootstrap

 Title   : move_id_to_bootstrap
 Usage   : $tree->move_id_to_bootstrap()
 Function: Move internal IDs to bootstrap slot
 Returns : undef
 Args    : undef

=cut

sub move_id_to_bootstrap{
    my $self = shift;
    foreach my $node ($self->root->nodes) {
	next if ($node->is_leaf);
	$node->bootstrap($node->id || '');
	$node->id('');
    }
}

=head2 to_file

 Title   : to_file
 Usage   : $tree->to_file()
 Function: Convenience method to write a tree to a file.
 Returns : nothing
 Args    : [Optional] A string indicating the format to use.

=cut
sub to_file {
    my $self = shift;
    my $file = shift;
    my $format = shift;

    $format = 'newick' unless (defined $format);
    my $iomod = "Bio::TreeIO::$format";
    $self->_load_module($iomod);

    my $io = $iomod->new(-format=>$format,-file=>">$file");
    $io->write_tree($self);
    $io->close;
}

1;
