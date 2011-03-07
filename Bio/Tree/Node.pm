#
# BioPerl module for Bio::Tree::Node
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Node - A Simple Tree Node

=head1 SYNOPSIS

    use Bio::Tree::Node;
    my $nodeL = Bio::Tree::Node->new(-id => 'L');
    my $nodeR = Bio::Tree::Node->new(-id => 'R');
    my $node = Bio::Tree::Node->new(-id => 'parent');
    $node->add_child($nodeL);
    $node->add_child($nodeR);

    print $node->ascii;

=head1 DESCRIPTION

The Node module is the stock class for representing nodes of a
phylogenetic tree in Bioperl.

The five NodeI methods (parent, children, add_child, remove_child,
branch_length) are implemented with parent / child links and branch
lengths stored within the object's hashtable. The order of child nodes
is maintained via nodes' internal_id values, allowing for the order of
child nodes to be flipped (e.g., node.reverse_subtree() ) or set
arbitrarily (e.g., node.set_child_order).

Metadata storage is enabled by inheriting the methods from
TagValueHolder. The node.bootstrap() and node.description() methods
just use pre-defined tags to store these common metadata types;
subclasses could add support for other common data types using this
approach.



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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Aaron Mackey, amackey-at-virginia-dot-edu
Sendu Bala,   bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::Node;
use vars qw($CREATIONORDER);
use strict;

use base ('Bio::Root::Root',
	  'Bio::Tree::NodeI',
	  'Bio::Tree::LabeledNodeI',
	  'Bio::Tree::NodeFunctionsI',
	  'Bio::Tree::TagValueHolder');

BEGIN {
    $CREATIONORDER = 1;
}


# Clones this node and the subtree below.
sub clone {
    my $self = shift;

    return _clone_and_add_children($self);
}

sub _clone_and_add_children {
    my $orig_node = shift;
    
    #my $clone_node = $orig_node->SUPER::clone;

    # Do a manual clone of all the node's object data.
    my $clone_node = new $orig_node;
    $clone_node->id($orig_node->id);
    $clone_node->branch_length($orig_node->branch_length);
    $clone_node->bootstrap($orig_node->bootstrap);
    my $tags = { %{$orig_node->{_tags}} };
    $clone_node->{_tags} = $tags;


    # Remove the cloned node from its parent (remove hanging ref)
    $clone_node->remove_from_parent;
    foreach my $clone_child ($clone_node->children) {
    	# Remove children from the cloned node (remove hanging refs)
    	$clone_node->remove_child($clone_child);
    }

    # Recurse, cloning children of the original node and adding them
    # to the current cloned node.
    foreach my $orig_child ($orig_node->children) {
	$clone_node->add_child(_clone_and_add_children($orig_child));
    }
    
    # Return the cloned node.
    return $clone_node;
}

=head2 new

 Title   : new
 Usage   : my $node = Bio::Tree::Node->new(-id => 'a',-branch_length => 0.1);
 Function: Builds a new Bio::Tree::Node object
 Returns : Bio::Tree::Node
 Args    : -id            => id / label for node [string] (optional but recommended)
           -descendents   => list of descendents [arrayref] (optional, they will be
                             updated so their ancestor pointer is this
                             node)
           -branch_length => branch length [integer] (optional)
           -bootstrap     => bootstrap value [string] (optional)
           -description   => description of the node [string] (optional)

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($children, $branchlen,$id,
      $bootstrap, $desc,$d) = $self->_rearrange([qw(
						    DESCENDENTS
						    BRANCH_LENGTH
						    ID
						    BOOTSTRAP
						    DESC
						    DESCRIPTION
						 )],
					     @args);
  $self->_register_for_cleanup(\&node_cleanup);
  $self->{'_desc'} = {}; # for descendents
  if( defined $d && defined $desc ) { 
      $self->warn("can only accept -desc or -description, not both, accepting -description");
      $desc = $d;
  } elsif( defined $d && ! defined $desc ) {
      $desc = $d;
  }
  defined $desc && $self->description($desc);
  defined $bootstrap && $self->bootstrap($bootstrap);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);
  if( defined $children ) {
      if( ref($children) !~ /ARRAY/i ) {
	  $self->warn("Must specify a valid ARRAY reference to initialize a Node's Descendents");
      }
      foreach my $c ( @$children ) { 	
	  $self->add_child($c);
      }
  }
  $self->sort_order($CREATIONORDER);
  $self->{'_internal_id'} = $CREATIONORDER++;
  return $self;
}


=head2 parent

 Title   : parent
 Usage   : my $parent = $node->parent;
 Function: Return this node's parent, or undef if this node has no
           parent. Implementation of the method defined in
           Bio::Tree::NodeI.
 Returns : Bio::Tree::NodeI or undef if no parent exists (i.e. the
           current node is the root)
 Args    : none

=cut

sub parent {
    my $self = shift;
    if (@_) {
        my $new_ancestor = shift;
        
        # we can set ancestor to undef
        if ($new_ancestor) {
            $self->throw("This is [$new_ancestor], not a Bio::Tree::NodeI")
		unless $new_ancestor->isa('Bio::Tree::NodeI');
        }
        
        my $old_ancestor = $self->{'_ancestor'} || '';
        if (!$old_ancestor || 
	    ($old_ancestor && ( !$new_ancestor || 
			       $new_ancestor ne $old_ancestor)) ) {
            if( $old_ancestor && ! $old_ancestor->{_removing_descendent}) {
		$old_ancestor->remove_Descendent($self);
	    }
            if ($new_ancestor && 
		! $new_ancestor->{_adding_descendent} ) { # avoid infinite recurse
                $self->{_setting_ancestor} = 1;
                $new_ancestor->add_Descendent($self, 1);
                $self->{_setting_ancestor} = 0;
            }
        }
        $self->{'_ancestor'} = $new_ancestor;
    }
    
    return $self->{'_ancestor'};
}


=head2 branch_length

 Title   : branch_length
 Usage   : $node->branch_length(0.5);
 Function: Get/Set the length of the branch between the current node
           and its parent. Implementation of the method defined in
           Bio::Tree::NodeI.
 Returns : Branch length
 Args    : new branch length (optional)

=cut

sub branch_length{
    my $self = shift;
    if( @_ ) {
	my $bl = shift;
	if( defined $bl &&
	    $bl =~ s/\[(\d+)\]// ) {
	    $self->bootstrap($1);
	}
	$self->{'_branch_length'} = $bl;
    }
    my $current_bl = $self->{'_branch_length'};
    #$current_bl = 1 unless (defined $current_bl && $current_bl =~ /^\-?\d+(\.\d+)?$/);
    return $current_bl;
}

=head2 children

 Title   : children
 Usage   : my @nodes = $node->children;
 Function: List the direct children of this node (but not all
           descendents; use .nodes() or .leaves() for this). Child
           nodes are sorted by their sort_order values. Implementation
           of the method defined in Bio::Tree::NodeI.
 Returns : Array of Bio::Tree::NodeI objects
 Args    : none

=cut

sub children {
   my ($self) = @_;

   return map { $_->[0] }
   sort { $a->[1] <=> $b->[1] } 
   map { [$_, $_->isa("Bio::Tree::Node") ? $_->sort_order : 0 ] }
   grep {defined $_}
   values %{$self->{'_desc'}};
}

=head2 add_child

 Title   : add_child
 Usage   : $node->add_child($new_node);
 Function: Add a child to a node. Implementation of the method defined
           in Bio::Tree::NodeI.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the node to add as a child of the current node

=cut

sub add_child{
   my ($self,$node,$ignoreoverwrite) = @_;
   return -1 if( ! defined $node );
   
   if( ! ref($node) ||
       ref($node) =~ /HASH/ ||
       ! $node->isa('Bio::Tree::NodeI') ) {
       $self->throw("Trying to add a Descendent who is not a Bio::Tree::NodeI");
       return -1;
   }
   
   $self->{_adding_descendent} = 1;
   # avoid infinite recurse
   $node->ancestor($self) unless $node->{_setting_ancestor}; 
   $self->{_adding_descendent} = 0;
   
   if( $self->{'_desc'}->{$node->internal_id} && ! $ignoreoverwrite ) {
       $self->throw("Going to overwrite a node which is $node that is already stored here, set the ignore overwrite flag (parameter 2) to true to ignore this in the future");
   }
   $self->{'_desc'}->{$node->internal_id} = $node;
   
}

=head2 remove_child

 Title   : remove_child
 Usage   : $node->remove_child($child_node);
 Function: Removes a child from the current node. Implementation of
           the method defined in Bio::Tree::NodeI.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the child to remove from the current node

=cut

sub remove_child{
   my ($self,@nodes) = @_;
   my $c= 0;
   foreach my $n ( @nodes ) { 
       if( $self->{'_desc'}->{$n->internal_id} ) {
        $self->{_removing_descendent} = 1;
        $n->ancestor(undef);
        $self->{_removing_descendent} = 0;
	   # should be redundant
	   $self->{'_desc'}->{$n->internal_id}->ancestor(undef);
	   delete $self->{'_desc'}->{$n->internal_id};
	   $c++;
       } else { 
	   if( $self->verbose ) {
	       $self->debug(sprintf("no node %s (%s) listed as a descendent in this node %s (%s)\n",$n->id, $n,$self->id,$self));
	       $self->debug("Descendents are " . join(',', keys %{$self->{'_desc'}})."\n");
	   }
       }
   }
   $c;
}

=head2 id

 Title   : id
 Usage   : $node->id($new_id_string);
 Function: Get/set the identifier / label for the node. Implementation
           of the method defined in Bio::Tree::LabeledNodeI.
 Returns : Value of the node's identifier
 Args    : String, new value for the node's ID

The id() method is designed to store and return any arbitrary
identifier string; thus, it does not produce properly Newick-escaped
labels. Use the L<id_output> method for creating Newick-compatible
strings.

=cut

sub id {
    my ($self, $value) = @_;
    if (defined $value) {
        $self->{'_id'} = $value;
    }
    my $val = $self->{'_id'};
    $val = '' if (!defined $val);
    return $val;
}

=head2 set_child_order

 Title   : set_child_order
 Usage   : $node->set_child_order($child_a,$child_b,$child_c);
 Function: Explicitly sets the sort order of the children beneath this node.
 Returns : nothing
 Args    : List of Bio::Tree::Node objects, which must contain all of the direct
           descendents of the current node

=cut

sub set_child_order {
    my $self = shift;
    my @nodes_in_new_order = @_;

    if (scalar(@nodes_in_new_order) != scalar($self->children)) {
	$self->warn("Number of nodes provided to set_child_order (".scalar(@nodes_in_new_order).") does not match the child count (".scalar($self->children).")! Doing nothing.");
	return;
    }

    # Collect the sort_orders of all the children nodes.
    my @sort_orders = map {$_->sort_order} $self->children;
    @sort_orders = sort {$a <=> $b}  @sort_orders;

    # Explicitly re-apply the sorted creation IDs to the new order of children
    foreach my $child (@nodes_in_new_order) {
	$child->sort_order(shift @sort_orders);
    }
}

=head2 reverse_children

 Title   : reverse_children
 Usage   : $node->reverse_children();
 Function: Reverse the order of the children directly beneath this node. Non-recursive.
 Returns : Nothing
 Args    : None

=cut

sub reverse_children {
    my $self = shift;

    my @children = $self->children;
    @children = reverse @children;
    $self->set_child_order(@children);
}

=head2 flip_subtree

 Title   : flip_subtree
 Usage   : $node->flip_subtree();
 Function: Flip the order of the entire sub-tree beneath this node (recursively).
 Returns : Nothing
 Args    : None

=cut

sub flip_subtree {
    my $self = shift;

    foreach my $node ($self->nodes) {
	$node->reverse_children;
    }
    $self->reverse_children;
}

=head2 internal_id

 Title   : internal_id
 Usage   : my $id = $node->internal_id;
 Function: Returns the unique internal ID of the current node. This ID
           is tied to Bioperl object creation and *will* change if the
           same tree is re-loaded from a file or string, but it is
           guaranteed to be unique for a given node in the tree.
 Returns : Integer ID
 Args    : None

=cut

sub internal_id {
    my $self = shift @_;
    return $self->{'_internal_id'} || 0;
}

=head2 sort_order

 Title   : sort_order
 Usage   : $node->sort_order(99);
 Function: Get/set the sort order for the current node. The sort order
           defines how a node's children are sorted when returned by
           the children() method. This method should not usually need
           to be called manually; instead, use methods like
           set_child_order or flip_subtree to handle common
           node-reordering operations.
 Returns : Integer, the sort order of the current node.
 Args    : Integer, the new sort order (optional)

=cut

sub sort_order {
    my $self = shift @_;
    $self->{'_sort_order'} = shift @_ if( @_);
    return $self->{'_sort_order'} || 0;
}

=head2 id_output

 Title   : id_output
 Usage   : my $nice_id = $node->id_output;
 Function: Return an id suitable for output in format like newick so
           that if it contains spaces or ():; characters it is
           properly quoted
 Returns : String
 Args    : none

From the Newick format description
(L<http://evolution.genetics.washington.edu/phylip/newicktree.html>):

"A name can be any string of printable characters except blanks,
colons, semicolons, parentheses, and square brackets. Because you may
want to include a blank in a name, it is assumed that an underscore
character ("_") stands for a blank; any of these in a name will be
converted to a blank when it is read in."

=cut

sub id_output{
    my $node = shift;
    my $id = $node->id;
    return unless( defined $id && length($id ) );
    # single quotes must become double quotes
    # $id =~ s/'/''/g;
    if( $id =~ /[\(\);:,\s]/ ) {
	$id = '"'.$id.'"';
    }
    return $id;
}

=head2 description

 Title   : description
 Usage   : $node->description("Yet Another node.");
 Function: Get / set the description for the current node.
 Returns : String
 Args    : String, the new description (optional)

=cut

sub description { 
    my $self = shift;
    if( @_ ) {
	if( $self->has_tag('description') ) {
	    $self->remove_tag('description');
	}
	$self->add_tag_value('description',shift);
    }
    return ($self->get_tag_values('description'))[0];
}

=head2 bootstrap

 Title   : bootstrap
 Usage   : $node->bootstrap(95);
 Function: Get / set the bootstrap support value for the current node
 Returns : Float
 Args    : Float, the new bootstrap value (optional)

=cut

sub bootstrap { 
    my $self = shift;
    if( @_ ) {
	if( $self->has_tag('B') ) {
	    $self->remove_tag('B');
	}
	$self->add_tag_value('B',shift);
    }
    return ($self->get_tag_values('B'))[0];
}

=head2 node_cleanup

 Title   : node_cleanup
 Usage   : $node->node_cleanup();
 Function: Cleans up this node by removing hanging references to child nodes.
 Returns : Nothing
 Args    : None

=cut

sub node_cleanup {
    my $self = shift;
    return unless defined $self;
    
    foreach my $child ($self->children) {
	$self->remove_child($child);
    }
    1;
}

1;
