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

    my $node = new Bio::Tree::Node();
    $node->add_Descendents($nodeL);
    $node->add_Descendents($nodeR);

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

use Bio::Root::Root;
use Bio::Tree::NodeI;

@ISA = qw(Bio::Root::Root Bio::Tree::NodeI);

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
  my ($children, $branchlen,$id,
      $bootstrap, $desc) = $self->_rearrange([qw(DESCENDENTS
						 BRANCH_LENGTH
						 ID
						 BOOTSTRAP
						 DESC
						 )],
					     @args);
  $self->{'_desc'} = {};
  defined $desc && $self->description($desc);
  defined $bootstrap && $self->bootstrap($bootstrap);
  defined $id && $self->id($id);
  defined $branchlen && $self->branch_length($branchlen);

  if( defined $children ) {
      if( ref($children) !~ /ARRAY/i ) {
	  $self->warn("Must specify a valid ARRAY reference to initialize a Node's Descendents");
      }
      foreach my $c ( @$children ) { 	
	  $self->add_Descendent($c);
      }
  }
  return $self;
}

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendant($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI

=cut

sub add_Descendent{
   my ($self,$node) = @_;
   return -1 if( ! defined $node ) ;
   if( ! $node->isa('Bio::Tree::NodeI') ) {
       $self->warn("Trying to add a Descendent who is not a Bio::Tree::NodeI");
       return -1;
   }
   # do we care about order?
   $node->ancestor($self);
   $self->{'_desc'}->{$node} = $node;
   return scalar keys %{$self->{'_desc'}};
}


=head2 each_Descendent

 Title   : each_Descendent
 Usage   : my @nodes = $node->each_Descendent;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : none

=cut

sub each_Descendent{
   my ($self) = @_;
   # order can be based on branch length (and sub branchlength)

   return sort { $a->height <=> $b->height } values %{$self->{'_desc'}};
}

=head2 get_Descendents

 Title   : get_Descendents
 Usage   : my @nodes = $node->get_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : none

=cut

=head2 ancestor

 Title   : ancestor
 Usage   : $obj->ancestor($newval)
 Function: Set the Ancestor
 Returns : value of ancestor
 Args    : newvalue (optional)

=cut

sub ancestor{
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

sub DESTROY {
    my ($self) = @_;
    # try to insure that everything is cleaned up
    $self->SUPER::DESTROY();
    if( defined $self->{'_desc'} &&
	ref($self->{'_desc'}) =~ /ARRAY/i ) {
	foreach my $n ( @{$self->{'_desc'}} ) {
	    $n->DESTROY();
	}
	$self->{'_desc'} = {};
    }
}

# The following methods are implemented by NodeI decorated interface

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf )
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf {
    my ($self) = @_;
    my $rc = 0;
    $rc = 1 if( ! defined $self->{'_desc'} ||	
		keys %{$self->{'_desc'}} == 0);
    return $rc;
}

=head2 to_string

 Title   : to_string
 Usage   : my $str = $node->to_string()
 Function: For debugging, provide a node as a string
 Returns : string
 Args    : none

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

1;

