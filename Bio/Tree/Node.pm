# $Id$
#
# BioPerl module for Bio::Tree::Node
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Node - A Simple Tree Node

=head1 SYNOPSIS

    my $node = new Bio::Tree::Node(-parent => $node,
				   -leaf   => 1);

    if( $node->is_leaf ) {
	# process a leafnode
    } else {
	# process a node
    }

=head1 DESCRIPTION

Makes a Tree Node suitable for building a node.

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::Node;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::Tree::NodeI;

@ISA = qw(Bio::Tree::NodeI Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Node();
 Function: Builds a new Bio::Tree::Node object 
 Returns : Bio::Tree::Node
 Args    : -parent => parent Node of this node
           -leaf   => boolean if this is a leaf node

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($parent,$leaf) = $self->_rearrange([qw(PARENT 
					     LEAF)],
					 @args);
  $self->set_parent($parent);
  $self->is_leaf($leaf);
  return $self;
}

=head2 add_child

 Title   : add_child
 Usage   : $node->add_child
 Function: Adds a node that is a child of a node,
           If this node already has a child will call add_child to its 
           child node
 Returns : integer - depth of the tree from this node to a leaf
 Args    : Bio::Tree::NodeI

=cut

sub add_child{
   my ($self,$node) = @_;   
   if( !defined $node || ! $node->isa('Bio::Tree::NodeI') ) {
       $self->throw("Must specify a Bio::Tree::NodeI to add_child");
   }
   my $r = 0;
   if( $self->is_leaf ) {
       $node->set_parent($self);
       $self->{'_child'} = $node;
   } else {
       $r = $self->{'_child'}->add_child($node);
   }
   return $r + 1;
}

=head2 get_child

 Title   : get_child
 Usage   : my $nodechild = $node->get_child
 Function: Retrieves the node child at this point
 Returns : Bio::Node::NodeI
 Args    : none


=cut

sub get_child{
   my ($self) = @_;
   return $self->{'_child'};
}

=head2 get_parent

 Title   : get_parent
 Usage   : my $node = $node->parent;
 Function: Gets a Node\'s parent node
 Returns : Null if this is top level node
 Args    : none

=cut

sub get_parent{
   my ($self,@args) = @_;
   return $self->{'_parent'};
}

=head2 set_parent

 Title   : set_parent
 Usage   : $obj->set_parent($newval)
 Function: Set the parent(s)
 Returns : value of set_parent
 Args    : newvalue (optional)

=cut

sub set_parent{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_parent'} = $value;
    }
    return $self->{'_parent'};
}

=head2 is_leaf

 Title   : is_leaf
 Usage   : if( $node->is_leaf ) 
 Function: Get/Set Leaf status
 Returns : boolean
 Args    : (optional) boolean

=cut

# implemented by NodeI decorated interface

1;
