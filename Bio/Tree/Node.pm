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
    
    use Bio::Tree::Node;
    my $nodeA = new Bio::Tree::Node();
    my $nodeL = new Bio::Tree::Node();
    my $nodeR = new Bio::Tree::Node();
    
    my $node = new Bio::Tree::Node(-left   => $nodeL,
				   -right  => $nodeR);

    $nodeA->set_Left_Descendent($node);
    
    print "node is not a leaf \n" if( $node->is_leaf);
   
=head1 DESCRIPTION

Makes a Tree Node suitable for building a Tree.

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
 Args    : -left   => pointer to Left descendent (optional)
           -right  => pointer to Right descenent (optional)
	   -branch_length => branch length [integer] (optional)
           -bootstrap => value   bootstrap value (string)
           -desc      => description of node
           -id        => unique id for node
=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($left,$right,$branchlen,$id,
      $bootstrap, $desc) = $self->_rearrange([qw(LEFT RIGHT
						 BRANCH_LENGTH
						 ID
						 BOOTSTRAP
						 DESC
						       )],
					     @args);

  defined $desc && $self->description($desc);  
  defined $bootstrap && $self->bootstrap($bootstrap);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);
  
  # these will operate fine if left or right is null
  $self->set_Left_Descendent($left);
  $self->set_Right_Descendent($right);
  return $self;
}

=head get_Left_Descendent

 Title   : get_Left_Descendent
 Usage   : my $left_D = $node->get_Left_Descendent 
 Function: Get the Left Descendant of a node
 Returns : Bio::Tree::NodeI
 Args    : none

=cut

sub get_Left_Descendent {
   my ($self) = @_;
   return $self->{'_left'};
}

=head get_Right_Descendent

 Title   : get_Right_Descendent
 Usage   : my $right_D = $node->get_Right_Descendent 
 Function: Get the Right Descendant of a node
 Returns : Bio::Tree::NodeI
 Args    : none

=cut

sub get_Right_Descendent {
   my ($self) = @_;
   return $self->{'_right'};
}

=head set_Left_Descendent

 Title   : set_Left_Descendent
 Usage   : $node->set_Left_Descendent($left) 
 Function: Set the Left Descendant of a node
 Returns : Bio::Tree::NodeI
 Args    : Bio::Tree::NodeI

=cut

sub set_Left_Descendent {
   my ($self,$value) = @_;
   if( defined $value ) {   
     if(! $value->isa('Bio::Tree::NodeI') ) { 
	 $self->throw("Must specify a valid Bio::Tree::NodeI when setting a descendent");
     } 
     
     $self->{'_left'} = $value;
     $self->{'_left'}->set_Ancestor($self);
   }
   return $self->{'_left'};
}

=head set_Right_Descendent

 Title   : set_Right_Descendent
 Usage   : $node->set_Right_Descendent($right) 
 Function: Set the Right Descendant of a node
 Returns : Bio::Tree::NodeI
 Args    : Bio::Tree::NodeI

=cut

sub set_Right_Descendent {
   my ($self,$value) = @_;
   if( defined $value ) {   
     if(! $value->isa('Bio::Tree::NodeI') ) { 
	 $self->throw("Must specify a valid Bio::Tree::NodeI when setting a descendent");
     } 
     $self->{'_right'} = $value;
     $self->{'_right'}->set_Ancestor($self);
   }
   return $self->{'_right'};
}

=head2 get_Ancestor

 Title   : get_Ancestor
 Usage   : my $node = $node->get_Ancestor;
 Function: Gets a Node\'s ancestor node
 Returns : Null if this is top level node
 Args    : none

=cut

sub get_Ancestor{
   my ($self,@args) = @_;
   return $self->{'_ancestor'};
}


=head2 set_Ancestor

 Title   : set_Ancestor
 Usage   : $obj->set_Ancestor($newval)
 Function: Set the Ancestor
 Returns : value of set_Ancestor
 Args    : newvalue (optional)

=cut

sub set_Ancestor{
   my ($self,$value) = @_;
   if( defined $value) {
       if(! $value->isa('Bio::Tree::NodeI') ) { 
	   $self->throw("Must specify a valid Bio::Tree::NodeI when setting the ancestor");
       }      
       $self->{'_ancestor'} = $value;
   }
   return $self->{'_ancestor'};
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
    if( defined $value) {
	$self->{'branch_length'} = $value;
    }
    return $self->{'branch_length'};
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
    if( defined $value ) {
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
   if( defined $value  ) {
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
   if( defined $value ) {
       $self->{'_id'} = $value;
   }
   return $self->{'_id'};
}

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf ) 
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

=head2 to_string

 Title   : to_string
 Usage   : my $str = $node->to_string()
 Function: For debugging, provide a node as a string
 Returns : string
 Args    : none


=cut

# implemented by NodeI decorated interface

1;
