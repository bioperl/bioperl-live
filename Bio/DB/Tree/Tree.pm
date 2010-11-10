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
    defined $id && $self->node_id($id);
    defined $label && $self->id($label);
    defined $root && $self->root_node($root);
    defined $store && $self->store($store);
    return $self;
}


=head2 tree_id

 Title   : tree_id
 Usage   : $obj->tree_id($newval)
 Function: DB related method signifying the internal creation order
 Returns : value of tree_id
 Args    : newvalue (optional)

=cut

sub tree_id {
    my $self = shift @_;
    $self->{'_tree_id'} = shift @_ if( @_);
    return $self->{'_tree_id'} || 0;
}

*_creation_id  = \&tree_id;

=head2 DB enabling methods

=head2 store

 Title   : store
 Usage   : $obj->store($newval)
 Function: Access to the L<Bio::DB::Tree::Store>
 Returns : value of node_id
 Args    : newvalue (optional)

=cut

sub store {
    my $self = shift;
    $self->{'_store'} = shift if( @_);
    return $self->{'_store'} || 0;
}

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
    $self->{'_id'} = shift;
    $self->_to_commit;
  } elsif ( ! exists $self->{'_label'} ) {
    $self->{'_label'} = $self->store->_fetch_node_label($self->node_id);
  }
  return $self->{'_label'};
}

=head2 root_node

 Title   : root_node
 Usage   : $obj->root_node($newval)
 Function: Get/set the Root Node for the tree (connection to the data)
 Returns : value of Bio::DB::Tree::Node
 Args    : newvalue (optional)

=cut

sub root_node {
  my $self = shift;
  if( @_ ) {
    $self->{'_rootnode'} = shift;
    $self->_to_commit;
  } elsif ( ! exists $self->{'_rootnode'} ) {
    $self->{'_rootnode'} = $self->store->_fetch_tree_root_node($self->tree_id);
  }
  return $self->{'_rootnode'};
}

=head2 root_node_id

 Title   : root_node_id
 Usage   : $obj->root_node_id
 Function: Get the root node id for the tree
 Returns : integer for the node
 Args    : None -- use root_node method to set the root node

=cut

sub root_node_id {
  shift->root_node->node_id;
}

sub _to_commit {
  my $self = shift;
  $self->{_commit_needed} = 1;
}


1;
