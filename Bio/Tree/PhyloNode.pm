# $Id$
#
# BioPerl module for Bio::Tree::PhyloNode
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::PhyloNode - A Phlyogenetic Tree Node 

=head1 SYNOPSIS

    my $phylonode = new Bio::Tree::PhyloNode(-parent => $node,
					     -leaf   => 1,
					     -branch_length => '0.301',
					     -boostrap=> '75', # 75% support 
					     -id     => 'ACT2_LYTPI',
					     -desc   => 'L.pictus Actin cytoskeleton 2');

    
=head1 DESCRIPTION

This object is a Phylogenetic Node suitable for building Phylogenetic trees.

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


package Bio::Tree::PhyloNode;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Tree::Node;

@ISA = qw(Bio::Tree::Node );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::PhyloNode();
 Function: Builds a new Bio::Tree::PhyloNode object 
 Returns : Bio::Tree::PhyloNode
 Args    : -parent => $parentNode (or null if root of tree)
           -leaf   => boolean    is a leaf node
           -bootstrap => value   bootstrap value (string)
           -desc      => description of node
           -id        => unique id for node

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($bootstrap, $desc, $id,
      $branchlen) = $self->_rearrange([qw(BOOTSTRAP
					  DESC
					  ID
					  BRANCH_LENGTH
					  )], @args);
  if( ! defined $id ) {
      $self->throw("Must define a valid ID for PhyloNode");
  }
  $desc && $self->description($desc);
  $bootstrap && $self->bootstrap($bootstrap);
  $id && $self->id($id);
  $branchlen && $self->branch_length($branchlen);
  return $self;

}

=head2 bootstrap

 Title   : bootstrap
 Usage   : $obj->bootstrap($newval)
 Function: 
 Example : 
 Returns : value of bootstrap
 Args    : newvalue (optional)


=cut

sub bootstrap{
   my ($self,$value) = @_;
   if( defined $value || ! defined $self->{'_boostrap'} ) {
       $value = '' unless defined $value;
       $self->{'_bootstrap'} = $value;
   }
   return $self->{'_bootstrap'};
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Example : 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description{
   my ($self,$value) = @_;
   if( defined $value || ! defined $self->{'_description'} ) {
       $value = '' unless defined $value;
       $self->{'_description'} = $value;
   }
   return $self->{'_description'};
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($self,$value) = @_;
   if( defined $value || ! defined $self->{'_id'} ) {
       $value = '' unless defined $value;
      $self->{'_id'} = $value;
    }
    return $self->{'_id'};

}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length($newval)
 Function: 
 Example : 
 Returns : value of branch_length
 Args    : newvalue (optional)


=cut

sub branch_length{
   my ($self,$value) = @_;
   if( defined $value || ! defined $self->{'_branch_length'} ) {
       $value = '' unless defined $value;
       $self->{'_branch_length'} = $value;
   }
   return $self->{'_branch_length'};
}

=head2 Bio::Tree::Node methods

=head2 is_leaf

 Title   : is_leaf
 Usage   : if( $node->is_leaf ) 
 Function: Get/Set Leaf status
           In Phylogenetic trees, leaf nodes are the only real nodes, the rest
           are inferred ancestors.
 Returns : boolean
 Args    : (optional) boolean

=head2 get_parent

 Title   : get_parent
 Usage   : my $node = $node->parent;
 Function: Gets a Node\'s parent node
 Returns : Null if this is top level node
 Args    : none

=head2 set_parent

 Title   : set_parent
 Usage   : $obj->set_parent($newval)
 Function: Set the parent(s)
 Returns : value of set_parent
 Args    : newvalue (optional)

=cut

1;
