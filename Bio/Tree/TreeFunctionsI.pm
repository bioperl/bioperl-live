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

Bio::Tree::TreeFunctionsI - Decorated Interface implementing basic Tree exploration methods

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = Bio::TreeIO->new(-format => 'newick', -file => 'tree.tre');

  my $tree = $in->next_tree;

  my @nodes = $tree->find_node('id1');

  if( $tree->is_monophyletic(-nodes => \@nodes, -outgroup => $outnode) ){
   #...
  }

=head1 DESCRIPTION

This interface provides a set of implementated Tree functions which
only use the defined methods in the TreeI or NodeI interface.

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

use base qw(Bio::Tree::TreeI);


sub print_tree { print shift->ascii }
sub ascii { shift->root->ascii }

sub get_lineage_nodes {
    my $self = shift;
    my $node = shift;
    return $node->lineage;
}

sub total_branch_length { shift->root->total_branch_length(@_) }
sub subtree_length { shift->root->subtree_length(@_) }
sub number_nodes { shift->root->node_count(@_) } # This alias sucks...

sub nodes {
    my $self = shift;

    # Slight difference between the Node 'nodes' method and Tree 'nodes' method --
    # the Tree version includes the root node in the returned array!

    return () unless ($self->root);
    my @nodes = $self->root->nodes;
    unshift @nodes, $self->root;
    return @nodes;
}
sub get_nodes { shift->nodes }
sub leaf_nodes { shift->root->leaves }
sub get_leaf_nodes { shift->root->leaves }

sub nodes_breadth_first { shift->root->nodes_breadth_first(@_) }
sub nodes_depth_first { shift->root->nodes_depth_first(@_) }

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
    # We must go root->leaf to get the correct answer to lca (in a world where
    # internal_id might not be uniquely assigned), but leaf->root is more
    # forgiving (eg. lineages may not all have the same root, or they may have
    # different numbers of 'minor' taxa inbeteen 'major' ones).
    #
    # I use root->leaf so that we can easily do multiple nodes at once - no
    # matter what taxa are below the lca, the lca and all its ancestors ought to
    # be identical.
    my @paths;
    foreach my $node (@nodes) {
	unless(ref($node) && $node->isa('Bio::Tree::NodeI')) {
	    $self->throw("Cannot process get_lca() with a non-NodeI object ($node)\n");
	}
        my @path = ($self->get_lineage_nodes($node), $node);
        push(@paths, \@path);
    }
    return unless @paths >= 2;
    my $lca;
    LEVEL: while ($paths[0] > 0) {
        my %node_ids;
        my $node;
        foreach my $path (@paths) {
            $node = shift(@{$path}) || last LEVEL;
            my $node_id = $node->internal_id;
            unless (defined $node_id) {
                $self->warn("One of the lineages had a node with no internal_id, can't calculate the common ancestor");
                return;
            }
            $node_ids{$node_id}++;
        }
        if (keys %node_ids == 1) {
            $lca = $node;
        }
        else {
            # at this point in the lineage the nodes are different; the previous
            # loop had the lca
            last LEVEL;
        }
    }
    # If the tree that we are contains the lca (get_lca could have been called
    # on an empty tree, since it works with plain Nodes), prefer to return the
    # node object that belongs to us
    if ($lca && $self->number_nodes > 0) {
        my $own_lca = $self->find_node(-internal_id => $lca->internal_id);
        $lca = $own_lca if $own_lca;
    }
    return $lca;
}

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
 Usage   : distance(-nodes => \@nodes )
 Function: returns the distance between two given nodes
 Returns : numerical distance
 Args    : -nodes => arrayref of nodes to test
           or ($node1, $node2)

=cut

sub distance {
    my ($self,@args) = @_;
    my ($nodes) = $self->_rearrange([qw(NODES)],@args);
    if( ! defined $nodes ) {
	$self->warn("Must supply two nodes or -nodes parameter to distance() method");
	return;
    }
    elsif (ref($nodes) eq 'ARRAY') {
	1;
    }
    elsif ( @args == 2) { # assume these are nodes...
	    $nodes = \@args;
    }
    else {
	$self->warn("Must supply two nodes or -nodes parameter to distance() method");
	return;
    }
    $self->throw("Must provide 2 nodes") unless @{$nodes} == 2;

    my $lca = $self->get_lca(@{$nodes});
    unless($lca) { 
        $self->warn("could not find the lca of supplied nodes; can't find distance either");
        return;
    }

    my $cumul_dist = 0;
    my $warned = 0;
    foreach my $current_node (@{$nodes}) {
        while (1) {
            last if $current_node eq $lca;
            if ($current_node->branch_length) {
                $cumul_dist += $current_node->branch_length;
            }
            elsif (! $warned) {
                $self->warn("At least some nodes do not have a branch length, the distance returned could be wrong");
                $warned = 1;
            }

            $current_node = $current_node->ancestor || last;
        }
    }

    return $cumul_dist;
}


sub is_monophyletic{
   my $self = shift;
   my $nodes = shift;
   my $outgroup = shift;

   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must supply nodes and outgroup parameters to the method
is_monophyletic");
       return;
   }
   if( ref($nodes) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid array reference for -nodes");
   }

   # Do a test of monophyly for the nodes specified
   # in comparison to a chosen outgroup

   my $clade_root = $self->get_lca(@{$nodes});
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return;
   }

   my $og_ancestor = $outgroup->ancestor;
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor->internal_id == $clade_root->internal_id ) {
           # monophyly is violated
           return 0;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }
   return 1;
}

sub is_paraphyletic{
    my $self = shift;
   my $nodes = shift;
   my $outgroup = shift;

   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must suply -nodes and -outgroup parameters to the method is_paraphyletic");
       return;
   }
   if( ref($nodes) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid array reference for -nodes");
       return;
   }

   # Tests whether or not a given set of nodes are paraphyletic
   # (representing the full clade) given an outgroup

   # Algorithm
   # Find the lca
   # Find all the nodes beneath the lca
   # Test to see that none are missing from the nodes list
   my %nodehash;
   foreach my $n ( @$nodes ) {
       $nodehash{$n->internal_id} = $n;
   }

   my $clade_root = $self->get_lca(@{$nodes});
   unless( defined $clade_root ) { 
       $self->warn("could not find clade root via lca");
       return;
   }

   my $og_ancestor = $outgroup->ancestor;

   # Is this necessary/correct for paraphyly test?
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor->internal_id == $clade_root->internal_id ) {
           # monophyly is violated, could be paraphyletic
           return -1;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }
   my $tree = Bio::Tree::Tree->new(-root     => $clade_root,
				  -nodelete => 1);

   foreach my $n ( $tree->get_nodes() ) { 
       next unless $n->is_Leaf();
       # if any leaf node is not in the list
       # then it is part of the clade and so the list
       # must be paraphyletic
       return 1 unless (  $nodehash{$n->internal_id} );
   }
   return 0;
}

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
    # For practical purposes you may want to use reroot_above instead.

    $new_root->parent->change_child_to_parent($new_root);
    $new_root->branch_length(undef);
    $self->set_root_node($new_root);

    return 1;
}

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
 Usage   : $tree->move_id_to_bootstrap
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

=head2 as_text

 Title   : as_text
 Usage   : my $tree_as_string = $tree->as_text($format)
 Function: Returns the tree as a string representation in the 
           desired format (currently 'newick', 'nhx', or 
           'tabtree')
 Returns : scalar string
 Args    : format type as specified by Bio::TreeIO
 Note    : This method loads the Bio::TreeIO::$format module
           on the fly, and commandeers the _write_tree_Helper
           routine therein to create the tree string. 

=cut

sub as_text {
    my $self = shift;
    my $format = shift;
    my $params_input = shift || {};

    $format = 'newick' unless (defined $format);

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


1;
