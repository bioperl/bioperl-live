#
# BioPerl module for Bio::DB::Tree::Node
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

Bio::DB::Tree::Node - A Database backed Tree Node

=head1 SYNOPSIS

    use Bio::DB::Tree::Node;
    my $nodeA = Bio::DB::Tree::Node->new();
    my $nodeL = Bio::DB::Tree::Node->new();
    my $nodeR = Bio::DB::Tree::Node->new();

    my $node = Bio::DB::Tree::Node->new();
    $node->add_Descendent($nodeL);
    $node->add_Descendent($nodeR);

    print "node is not a leaf \n" if( $node->is_leaf);

=head1 DESCRIPTION

Makes a Tree Node suitable for building a Tree.

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Tree::Node;
use strict;

use base qw(Bio::Root::Root Bio::Tree::NodeI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Node->new();
 Function: Builds a new Bio::Tree::Node object
 Returns : Bio::DB::Tree::Node
 Args    : -node_id            => internal id
           -parent_id          => parent id
           -id                 => Node label ID
           -branch_length      => branch length
           -store              => L<Bio::DB::Tree::Store> object

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($id,$label,$bl,$parent,$store) = $self->_rearrange([qw(NODE_ID
							   ID
							   BRANCH_LENGTH
							   PARENT_ID
							   STORE
							 )],
						       @args);
    defined $id && $self->node_id($id);
    defined $label && $self->id($label);
    defined $bl && $self->branch_length($bl);
    defined $parent && $self->parent_id($parent);
    defined $store && $self->store($store);
    return $self;
}


=head2 node_id

 Title   : node_id
 Usage   : $obj->node_id($newval)
 Function: DB related method signifying the internal creation order
 Returns : value of node_id
 Args    : newvalue (optional)

=cut

sub node_id {
    my $self = shift;
    return $self->{'_node_id'} = shift @_ if( @_);
    return $self->{'_node_id'};
}

=head2 parent_id

 Title   : parent_id
 Usage   : $obj->parent_id($newval)
 Function: DB related method signifying the parent id
 Returns : value of parent_id
 Args    : newvalue (optional)

=cut

sub parent_id {
    my $self = shift;
    if ( @_ ) {
        $self->{'_parent_id'} = shift @_;
        delete $self->{'_ancestor'};
        $self->_dirty(1);
    } elsif ( ! exists($self->{'_parent_id'})) {
        $self->_load_from_db;
    }
    return $self->{'_parent_id'};
}

*_creation_id  = \&node_id;

=head2 DB enabling methods

=head2 store

 Title   : store
 Usage   : $obj->store($newval)
 Function: Access to the L<Bio::DB::Tree::Store>
 Returns : the database store to which this object is connected
 Args    : newvalue (optional)

=cut

sub store {
    my $self = shift;
    $self->{'_store'} = shift if( @_);
    return $self->{'_store'} || 0;
}

=head2 save

 Title   : save
 Usage   :
 Function: Save the current state of the object to the database (if it
           has diverged).
 Example :
 Returns : True on success and false otherwise.
 Args    : none

=cut

sub save{
    my $self = shift;
    my $rv = 1;
    if ($self->_dirty) {
        $rv = $self->store->update_node($self);
        $self->_dirty(0) if $rv;
    }
    return $rv;
}

=head2 Bio::Tree::NodeI methods

=head2 create_node_on_branch

 Title   : create_node_on_branch
 Usage   : $node->create_node_on_branch($at_length)
 Function: Create a node on the ancestral branch of the calling
           object. 
 Example :
 Returns : the created node
 Args    : -POSITION=>$absolute_branch_length_from_caller (default)
           -FRACTION=>$fraction_of_branch_length_from_caller
           -ANNOT=>{ -id => "the id", -desc => "the description" }
           -FORCE, set to allow nodes with zero branch lengths

=cut

sub create_node_on_branch{
   my ($self,@args) = @_;
   my ($pos, $frac, $annot, $force) = $self->_rearrange([qw(POSITION 
							    FRACTION ANNOT 
							    FORCE)], @args);
   # are we creating this?
}

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendent($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.

=cut

sub add_Descendent{
  # are we doing this ?
  shift->throw("not implemented");
}

=head2 each_Descendent

 Title   : each_Descendent($sortby)
 Usage   : my @nodes = $node->each_Descendent;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : $sortby [optional] "height", "creation", "alpha", "revalpha",
           or coderef to be used to sort the order of children nodes.

=cut

sub each_Descendent{
  my $self = shift;
  my $sortby = shift;
  $self->store->_fetch_node_children($self->node_id, $sortby);
}


=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $node->remove_Descendent($node_foo);
 Function: Removes a specific node from being a Descendent of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects which have been previously
           passed to the add_Descendent call of this object.

=cut

sub remove_Descendent{
  my $self = shift;
  $self->store->_remove_node($self->node_id);
}


=head2 remove_all_Descendents

 Title   : remove_all_Descendents
 Usage   : $node->remove_All_Descendents()
 Function: Cleanup the node's reference to descendents and reset
           their ancestor pointers to undef, if you don't have a reference
           to these objects after this call they will be cleaned up - so
           a get_nodes from the Tree object would be a safe thing to do first
 Returns : nothing
 Args    : none

=cut

sub remove_all_Descendents{
  my $self = shift;
    $self->store->_remove_node($self->node_id);

}

=head2 get_all_Descendents

 Title   : get_all_Descendents
 Usage   : my @nodes = $node->get_all_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : none

=cut

sub get_all_Descendents {
  my $self = shift;
  my $sortby = shift;
  $self->store->_fetch_node_all_children($self->node_id, $sortby);
}

=head2 ancestor

 Title   : ancestor
 Usage   : $obj->ancestor($newval)
 Function: Get or set the Ancestor (parent node)
 Returns : Parent node
 Args    : on set, new value (optional)

=cut

sub ancestor {
    my $self = shift;

    if ( @_ ) {
        my $parent = shift;
        if (ref($parent)) {
            $self->{'_ancestor'} = $parent;
        } else {
            # if it's not a reference, assume it's a primary key
            $self->parent_id($parent);
            $self->{'_ancestor'} = $self->store->fetch_node($parent);
        }
        $self->_dirty(1);
    } elsif (! exists($self->{'_ancestor'})) {
        $self->{'_ancestor'} = $self->store->fetch_node($self->parent_id);
    }
    return $self->{'_ancestor'};
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: Get/Set the branch length
 Returns : value of branch_length
 Args    : newvalue (optional)

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
    $self->invalidate_height();
    $self->_dirty(1);
  } elsif ( ! exists $self->{'_branch_length'} ) {
      $self->_load_from_db;
  }
  return $self->{'_branch_length'};
}

=head2 bootstrap

 Title   : bootstrap
 Usage   : $obj->bootstrap($newval)
 Function: Get/Set the bootstrap value
 Returns : value of bootstrap
 Args    : newvalue (optional)

=cut

sub bootstrap {

}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: Get/Set the description string
 Returns : value of description
 Args    : newvalue (optional)

=cut

sub description {
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: The human readable identifier for the node 
 Returns : value of human readable id
 Args    : newvalue (optional)

"A name can be any string of printable characters except blanks,
colons, semicolons, parentheses, and square brackets. Because you may
want to include a blank in a name, it is assumed that an underscore
character ("_") stands for a blank; any of these in a name will be
converted to a blank when it is read in."  

from L<http://evolution.genetics.washington.edu/phylip/newicktree.html>

Also note that these objects now support spaces, ();: because we can
automatically quote the strings if they contain these characters.  The
L<id_output> method does this for you so use the id() method to get
the raw string while L<id_output> to get the pre-escaped string.

=cut

sub id {
  my $self = shift;
  if( @_ ) {
    $self->{'_label'} = shift;
    $self->_dirty(1);
  } elsif ( ! exists $self->{'_label'} ) {
      $self->_load_from_db;      
  }
  return $self->{'_label'};
}

=head2 Helper Functions

=cut

=head2 id_output

 Title   : id_output
 Usage   : my $id = $node->id_output;
 Function: Return an id suitable for output in format like newick
           so that if it contains spaces or ():; characters it is properly 
           quoted
 Returns : $id string if $node->id has a value
 Args    : none

=cut

# implemented in NodeI interface 


=head2 internal_id

 Title   : internal_id
 Usage   : my $internalid = $node->internal_id
 Function: Returns the internal unique id for this Node
           (a monotonically increasing number for this in-memory implementation
            but could be a database determined unique id in other 
	    implementations)
 Returns : unique id
 Args    : none

=cut

sub internal_id {
    shift->node_id(@_);
}

=head2 Bio::Node::NodeI decorated interface implemented

The following methods are implemented by L<Bio::Node::NodeI> decorated
interface.

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf )
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf {
  my $self = shift;
  return $self->store->_is_Leaf($self->node_id);
}

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength to get to the tip.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

sub height {
}

=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the node height in the tree
 Returns : nothing
 Args    : none

=cut

sub invalidate_height { 
}

=head2 set_tag_value

 Title   : set_tag_value
 Usage   : $node->set_tag_value($tag,$value)
           $node->set_tag_value($tag,@values)
 Function: Sets a tag value(s) to a node. Replaces old values.
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub set_tag_value{

}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $node->add_tag_value($tag,$value)
 Function: Adds a tag value to a node 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub add_tag_value{
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $node->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag {
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $node->remove_all_tags()
 Function: Removes all tags
 Returns : None
 Args    : None

=cut

sub remove_all_tags{
}


=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $node->get_all_tags()
 Function: Gets all the tag names for this Node
 Returns : Array of tagnames
 Args    : None

=cut

sub get_all_tags{
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $node->get_tag_values($tag)
 Function: Gets the values for given tag ($tag)
 Returns : In array context returns an array of values
           or an empty list if tag does not exist.
           In scalar context returns the first value or undef.
 Args    : $tag - tag name

=cut

sub get_tag_values{
}

=head2 has_tag

 Title   : has_tag
 Usage   : $node->has_tag($tag)
 Function: Boolean test if tag exists in the Node
 Returns : Boolean
 Args    : $tag - tagname

=cut

sub has_tag {
}


=head2 reverse_edge

 Title   : reverse_edge
 Usage   : $node->reverse_edge(child);
 Function: makes child be a parent of node
 Requires: child must be a direct descendent of node
 Returns : 1 on success, 0 on failure
 Args    : Bio::Tree::NodeI that is in the tree

=cut

sub reverse_edge {
}

=head1 Private methods

=head2 _dirty

 Title   : _dirty
 Usage   : $obj->_dirty($newval)
 Function: Whether or not the object is "dirty", i.e., may have a
           state that is diverged from the state in the corresponding
           database record, and to which therefore the database needs
           to be updated.
 Example : 
 Returns : True if object is dirty and false otherwise.
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _dirty{
    my $self = shift;

    return $self->{'_dirty'} = shift if @_;
    return $self->{'_dirty'};
}

=head2 _load_from_db

 Title   : _load_from_db
 Usage   : $obj->_load_from_db
 Function: Loads the object's state from the database if it hasn't
           been loaded before. Even after the state has been loaded,
           there may be additional lazy-loaded attributes that are
           loaded separately only once they are needed.
 Example : 
 Returns : True on success and false otherwise.
 Args    : none

=cut

sub _load_from_db{
    my $self = shift;
    my $rv = 1;
    if (! $self->{'_loaded'}) {
        $rv = $self->store->populate_node($self);
        $self->{'_loaded'} = 1 if $rv;
    }
    return $rv;
}

1;

