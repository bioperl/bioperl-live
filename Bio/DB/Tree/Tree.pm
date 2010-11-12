#
# BioPerl module for Bio::DB::Tree::Tree
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

Bio::DB::Tree::Tree - A Database backed Tree object

=head1 SYNOPSIS

    use Bio::DB::Tree::Tree;
    use Bio::DB::Tree::Node;
    my $nodeA = Bio::DB::Tree::Node->new();
    my $nodeL = Bio::DB::Tree::Node->new();
    my $nodeR = Bio::DB::Tree::Node->new();

    my $node = Bio::DB::Tree::Node->new();
    $node->add_Descendent($nodeL);
    $node->add_Descendent($nodeR);

    my $tree = Bio::DB::Tree::Tree->new(-root_node => $nodeA);

=head1 DESCRIPTION

Makes a Tree object suitable for processing DB backed tree data.

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

package Bio::DB::Tree::Tree;
use strict;

use base qw(Bio::Root::Root Bio::Tree::TreeI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Tree->new();
 Function: Builds a new Bio::Tree::Tree object
 Returns : Bio::DB::Tree::Tree
 Args    : -tree_id            => internal database id
           -id                 => tree name
           -root_node          => root node for tree
           -store              => L<Bio::DB::Tree::Store> object

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($id,$label,$root,$store) = $self->_rearrange([qw(TREE_ID
						   ID
						   ROOT_NODE
						   STORE
						 )],
					       @args);
    defined $id && $self->tree_id($id);
    defined $label && $self->id($label);
    defined $root && $self->root($root);
    defined $store && $self->store($store);
    $self->_dirty(0);
    return $self;
}

=head2 DB enabling methods

=head2 tree_id

 Title   : tree_id
 Usage   : $obj->tree_id($newval)
 Function: The primary key of the tree in the underlying storage.
 Returns : value of tree_id
 Args    : newvalue (optional)

=cut

sub tree_id {
    my $self = shift;
    return $self->{'_tree_id'} = shift if @_;
    return $self->{'_tree_id'};
}

=head2 store

 Title   : store
 Usage   : $obj->store($newval)
 Function: Access to the L<Bio::DB::Tree::Store>
 Returns : the database store object
 Args    : on set, newvalue (optional)

=cut

sub store {
    my $self = shift;
    return $self->{'_store'} = shift if @_;
    return $self->{'_store'};
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
        $rv = $self->store->update_tree($self);
        $self->_dirty(0) if $rv;
    }
    return $rv;
}

=head2 Tree properties

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: A string identifier for the tree
 Returns : value of readable id/accession for a tree
 Args    : newvalue (optional)

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

=head2 root

 Title   : root
 Usage   : $obj->root($newval)
 Function: Get/set the Root Node for the tree (connection to the data)
 Returns : A Bio::Tree::NodeI object (here, implemented by Bio::DB::Tree::Node)
 Args    : newvalue (optional)

=cut

sub root {
  my $self = shift;
  if( @_ ) {
    $self->{'_rootnode'} = shift;
    $self->_dirty(1);
  } elsif ( ! exists $self->{'_rootnode'} ) {
    $self->_load_from_db;
  }
  return $self->{'_rootnode'};
}

=head2 rooted

 Title   : rooted
 Usage   : $tree->rooted(0);
 Function: Get/Set Boolean whether or not this tree should be considered
           to be rooted or not. Note that the TreeI method is_rooted
           is an interface to this getter.
 Returns : boolean
 Args    : on set, new boolean value

=cut

sub rooted {
    my $self = shift;
    if( @_ ) {
        $self->{'_rooted'} = shift;
        $self->_dirty(1);
    }
    return $self->{'_rooted'};
}

=head2 is_rooted

 Title   : is_rooted
 Usage   : die unless ($tree->is_rooted);
 Function: Indicates whether this should be handled as a rooted tree.
 Returns : 1 if this tree is rooted; 0 if unrooted.
 Args    : none

=cut

sub is_rooted {
    my $self = shift;
    if ( ! exists $self->{'_rooted'} ) {
        $self->_load_from_db;
    }
    return $self->{'_rooted'} if defined $self->{'_rooted'};
    # Default to a rooted tree.
    return 1;	
}

=head2 score

 Title   : score
 Usage   : $score = $obj->score()
 Function: Get (or set) the score of the tree.
 Example : 
 Returns : value of score (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub score{
    my $self = shift;

    if( @_ ) {
        $self->{'_score'} = shift;
        $self->_dirty(1);
    } elsif ( ! exists $self->{'_score'} ) {
        $self->_load_from_db;
    }
    return $self->{'_score'};
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
        $self->{'_loaded'} = 1; # stop recursion while we're doing this
        $rv = $self->store->populate_tree($self);
        $self->{'_loaded'} = $rv;
    }
    return $rv;
}

1;
