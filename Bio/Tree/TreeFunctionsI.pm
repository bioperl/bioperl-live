# $Id$
#
# BioPerl module for Bio::Tree::TreeFunctionsI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::TreeFunctionsI - Decorated Interface implementing basic Tree exploration methods

=head1 SYNOPSIS

  use Bio::TreeIO;
  my $in = new Bio::TreeIO(-format => 'newick', -file => 'tree.tre');

  my $tree = $in->next_tree;

  my @nodes = $tree->find_nodes('id1');

  if( $tree->is_monophyletic(-clade => @nodes, -outgroup => $outnode) ){

  }

=head1 DESCRIPTION

Describe the interface here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich, Aaron Mackey, Justin Reese

Email jason@bioperl.org
Email amackey@virginia.edu
Email jtr4v@virginia.edu

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreeFunctionsI;
use vars qw(@ISA);
use strict;
use  Bio::Tree::TreeI;

@ISA = qw(Bio::Tree::TreeI);

=head2 find_node

 Title   : find_node
 Usage   : my @nodes = $self->find_node(-id => 'node1');
 Function: returns all nodes that match a specific field, by default this
           is id, but different branch_length, 
 Returns : List of nodes which matched search
 Args    : text string to search for
           OR
           -fieldname => $textstring

=cut

sub find_node {
   my ($self,$type,$field) = @_;
   if( ! defined $type ) { 
       $self->warn("Must request a either a string or field and string when searching");
   }

   # all this work for a '-' named field
   # is so that we could potentially 
   # expand to other constraints in 
   # different implementations
   # like 'find all nodes with boostrap < XX'

   if( ! defined $field ) { 
       # only 1 argument, default to searching by id
       $field= $type; 
       $type = 'id';
   } else {   
       $type =~ s/^-//;
   }

   # could actually do this by testing $rootnode->can($type) but
   # it is possible that a tree is implemeted with different node types
   # - although it is unlikely that the root node would be richer than the 
   # leaf nodes.  Can't handle NHX tags right now

   unless( $type eq 'id' || $type eq 'name' ||
	   $type eq 'bootstrap' || $type eq 'description' ||
	   $type eq 'internal_id') {
       $self->warn("unknown search type $type - will try anyways");
   } 
   my @nodes = grep { $_->can($type) && defined $_->$type() &&
		     $_->$type() eq $field } $self->get_nodes();
   
   if ( wantarray) { 
       return @nodes;
   } else { 
       if( @nodes > 1 ) { 
	   $self->warn("More than 1 node found but caller requested scalar, only returning first node");
       }
       return shift @nodes;
   }
}

=head2 remove_Node

 Title   : remove_Node
 Usage   : $tree->remove_Node($node)
 Function: Removes a node from the tree
 Returns : boolean represent status of success
 Args    : either Bio::Tree::NodeI or string of the node id


=cut

sub remove_Node {
   my ($self,$input) = @_;
   my $node = undef;
   unless( ref($input) ) {
       $node = $self->find_node($input);
   }  elsif( ! $input->isa('Bio::Tree::NodeI') ) {
       $self->warn("Did not provide either a valid Bio::Tree::NodeI object to remove_node or the node name");
       return 0;
   } else { 
       $node = $input;
   }
   if( ! $node->ancestor && $self->get_root_node->internal_id != $node->internal_id) {
       $self->warn("Node (".$node->to_string . ") has no ancestor, can't remove!");
   } else { 
       $node->ancestor->remove_Descendent($node);
   }
}


# Added for Justin Reese by Jason

=head2 get_lca

 Title   : get_lca
 Usage   : get_lca(-nodes => \@nodes )
 Function: given two nodes, returns the lowest common ancestor
 Returns : node object
 Args    : -nodes => arrayref of nodes to test


=cut

sub get_lca {
    my ($self,@args) = @_;
    my ($nodes) = $self->_rearrange([qw(NODES)],@args);
   if( ! defined $nodes ) {
       $self->warn("Must supply -nodes parameter to get_lca() method");
       return undef;
   }
    my ($node1,$node2) = $self->_check_two_nodes($nodes);
    return undef unless $node1 && $node2;

    # algorithm: Start with first node, find and save every node from it to
    #    root. Then start with second node; for it and each of its ancestor
    #    nodes, check to see if it's in the first node's ancestor list - if
    #    so it is the lca.
    #
    # This is very slow and naive, but I somehow doubt the overhead
    # of mapping the tree to a complete binary tree and doing the linear
    # lca search would be worth the overhead, especially for small trees.
    # Maybe someday I'll write a linear get_lca and find out.

    # find and save every ancestor of node1 (including itself)

    my %node1_ancestors;	# keys are internal ids, values are objects
    my $place = $node1;		# start at node1

    while ( $place ){
	$node1_ancestors{$place->internal_id} = $place;
	$place = $place->ancestor;
    }

    # now climb up node2, for each node checking whether 
    # it's in node1_ancestors
    $place = $node2;		# start at node2
    while ( $place ){
	foreach my $key ( keys %node1_ancestors ){ # ugh
	    if ( $place->internal_id == $key){
		return $node1_ancestors{$key};
	    }
	}
	$place = $place->ancestor;
    }
    $self->warn("Could not find lca!"); # should never execute, 
                                        # if so, there's a problem
    return undef;
}

# Added for Justin Reese by Jason

=head2 distance

 Title   : distance
 Usage   : distance(-nodes => \@nodes )
 Function: returns the distance between two given nodes
 Returns : numerical distance
 Args    : -nodes => arrayref of nodes to test


=cut

sub distance {
    my ($self,@args) = @_;
    my ($nodes) = $self->_rearrange([qw(NODES)],@args);
    if( ! defined $nodes ) {
	$self->warn("Must supply -nodes parameter to distance() method");
	return undef;
    }
    my ($node1,$node2) = $self->_check_two_nodes($nodes);
    # algorithm:

    # Find lca: Start with first node, find and save every node from it
    # to root, saving cumulative distance. Then start with second node;
    # for it and each of its ancestor nodes, check to see if it's in
    # the first node's ancestor list - if so it is the lca. Return sum
    # of (cumul. distance from node1 to lca) and (cumul. distance from
    # node2 to lca)

    # find and save every ancestor of node1 (including itself)

    my %node1_ancestors;	# keys are internal ids, values are objects
    my %node1_cumul_dist;	# keys are internal ids, values 
    # are cumulative distance from node1 to given node
    my $place = $node1;		# start at node1
    my $cumul_dist = 0;

    while ( $place ){
	$node1_ancestors{$place->internal_id} = $place;
	$node1_cumul_dist{$place->internal_id} = $cumul_dist;
	if ($place->branch_length) {
	    $cumul_dist += $place->branch_length; # include current branch
	                                          # length in next iteration
	}
	$place = $place->ancestor;
    }

    # now climb up node2, for each node checking whether 
    # it's in node1_ancestors
    $place = $node2;  # start at node2
    $cumul_dist = 0;
    while ( $place ){
	foreach my $key ( keys %node1_ancestors ){ # ugh
	    if ( $place->internal_id == $key){ # we're at lca
		return $node1_cumul_dist{$key} + $cumul_dist;
	    }
	}
	# include current branch length in next iteration
	$cumul_dist += $place->branch_length; 
	$place = $place->ancestor;
    }
    $self->warn("Could not find distance!"); # should never execute, 
    # if so, there's a problem
    return undef;
}

# helper function to check lca and distance arguments

sub _check_two_nodes {    
    my ($self, $nodes) = @_;

   if( ref($nodes) !~ /ARRAY/i ||
       !ref($nodes->[0]) ||
       !ref($nodes->[1])
       ) {
       $self->warn("Must provide a valid array reference for -nodes");
       return undef;
   } elsif( scalar(@$nodes) > 2 ){
       $self->warn("More than two nodes given, using first two");
   } elsif( scalar(@$nodes) < 2 ){
       $self->warn("-nodes parameter does not contain reference to two nodes");
       return undef;
   }    
    unless( $nodes->[0]->isa('Bio::Tree::NodeI') &&
	    $nodes->[1]->isa('Bio::Tree::NodeI') ) {
	$self->warn("Did not provide valid Bio::Tree::NodeI objects as nodes\n");
	return undef;
    }
    return @$nodes;
}


=head2 is_monophyletic

 Title   : is_monophyletic
 Usage   : if( $tree->is_monophyletic(-nodes => \@nodes, 
				      -outgroup => $outgroup)
 Function: Will do a test of monophyly for the nodes specified
           in comparison to a chosen outgroup
 Returns : boolean
 Args    : -nodes => arrayref of nodes to test
           -outgroup => outgroup to serve as a reference


=cut

sub is_monophyletic{
   my ($self,@args) = @_;
   my ($nodes,$outgroup) = $self->_rearrange([qw(NODES OUTGROUP)],@args);

   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must supply -nodes and -outgroup parameters to the method
is_monophyletic");
       return undef;
   }
   if( ref($nodes) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid array reference for -nodes");
   }
   my $clade_root;
   # this is to combine multiple tests into a single node
   # order doesn't really matter as long as get_lca does its job right
   while( @$nodes > 2 ) { 
       my ($a,$b) = ( shift @$nodes, shift @$nodes);
       $clade_root = $self->get_lca(-nodes => [$a,$b] );
       unshift @$nodes, $clade_root;
   }
   $clade_root = $self->get_lca(-nodes => $nodes );
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

=head2 is_paraphyletic

 Title   : is_paraphyletic
 Usage   : if( $tree->is_paraphyletic(-nodes =>\@nodes,
				      -outgroup => $node) ){ }
 Function: Tests whether or not a given set of nodes are paraphyletic
           (representing the full clade) given an outgroup
 Returns : [-1,0,1] , -1 if the group is not monophyletic
                       0 if the group is not paraphyletic
                       1 if the group is paraphyletic
 Args    : -nodes => Array of Bio::Tree::NodeI objects which are in the tree
           -outgroup => a Bio::Tree::NodeI to compare the nodes to


=cut

sub is_paraphyletic{
   my ($self,@args) = @_;
   my ($nodes,$outgroup) = $self->_rearrange([qw(NODES OUTGROUP)],@args);
   
   if( ! defined $nodes || ! defined $outgroup ) {
       $self->warn("Must suply -nodes and -outgroup parameters to the method is_paraphyletic");
       return undef;
   }
   if( ref($nodes) !~ /ARRAY/i ) { 
       $self->warn("Must provide a valid array reference for -nodes");
   }

   # Algorithm
   # Find the lca
   # Find all the nodes beneath the lca
   # Test to see that none are missing from the nodes list
   my %nodehash;
   foreach my $n ( @$nodes ) {
       $nodehash{$n->internal_id} = $n;
   }
   while( @$nodes > 2 ) { 
       unshift @$nodes, $self->get_lca(-nodes => [( shift @$nodes, 
						    shift @$nodes)] );
   }
   my $clade_root = $self->get_lca(-nodes => $nodes );
   my $og_ancestor = $outgroup->ancestor;

   # Is this necessary/correct for paraphyly test?
   while( defined ($og_ancestor ) ) {
       if( $og_ancestor->internal_id == $clade_root->internal_id ) {
           # monophyly is violated, could be paraphyletic
           return -1;
       }
       $og_ancestor = $og_ancestor->ancestor;
   }
   my $tree = new Bio::Tree::Tree(-root => $clade_root);
   
   foreach my $n ( $tree->get_nodes() ) { 
       next unless $n->is_Leaf();
       # if any leaf node is not in the list
       # then it is part of the clade and so the list
       # must be paraphyletic
       return 1 unless (  $nodehash{$n->internal_id} );
   }
   return 0;
}


=head2 reroot

 Title   : reroot_tree
 Usage   : $tree->reroot($node);
 Function: Reroots a tree either making a new node the root
 Returns : 1 on success, 0 on failure
 Args    : Bio::Tree::NodeI that is in the tree

o
=cut

sub reroot {
    my ($self,$node) = @_;
    unless (defined $node && $node->isa("Bio::Tree::NodeI")) {
	$self->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
	return 0;
    }
    # Still need to validate that $node is actually in the tree....
    if( $node->is_Leaf() ) {
	$self->warn("Asking to root with a leaf, will use the leaf's ancestor");
	$node = $node->ancestor;
    }

    my $current_root = $self->get_root_node;
    my $old_ancestor = $node->ancestor;

    if( ! $old_ancestor ) {
	$self->warn("Node requested for reroot is already the root node!");
	return 0;
    }

    $old_ancestor->remove_Descendent($node);
    $node->add_Descendent($current_root);
    $self->set_root_node($node);
    return 1;
}

1;
