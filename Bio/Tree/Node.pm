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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Aaron Mackey amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tree::Node;
use vars qw(@ISA $CREATIONORDER);
use strict;

use Bio::Root::Root;
use Bio::Tree::NodeI;

@ISA = qw(Bio::Root::Root Bio::Tree::NodeI);

BEGIN { 
    $CREATIONORDER = 0;
}


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Node();
 Function: Builds a new Bio::Tree::Node object
 Returns : Bio::Tree::Node
 Args    : -left          => pointer to Left descendent (optional)
           -right         => pointer to Right descenent (optional)
	   -branch_length => branch length [integer] (optional)
           -bootstrap     => value   bootstrap value (string)
           -description   => description of node
           -id            => unique id for node

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($children, $branchlen,$id,
      $bootstrap, $desc,$d) = $self->_rearrange([qw(DESCENDENTS
						 BRANCH_LENGTH
						 ID
						 BOOTSTRAP
						 DESC
						 DESCRIPTION
						 )],
					     @args);
  $self->{'_desc'} = {};
  if( $d && $desc ) { 
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
	  $self->add_Descendent($c);
      }
  }
  $self->_creation_id($CREATIONORDER++);
  return $self;
}

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendant($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.

=cut

sub add_Descendent{
   my ($self,$node,$ignoreoverwrite) = @_;
   return -1 if( ! defined $node ) ;
   if( ! $node->isa('Bio::Tree::NodeI') ) {
       $self->warn("Trying to add a Descendent who is not a Bio::Tree::NodeI");
       return -1;
   }
   # do we care about order?
   $node->ancestor($self);
   if( $self->{'_desc'}->{$node->internal_id} && ! $ignoreoverwrite ) {
       $self->throw("Going to overwrite a node which is $node that is already stored here, set the ignore overwrite flag (parameter 2) to true to ignore this in the future");
   }
   
   $self->{'_desc'}->{$node->internal_id} = $node; # is this safely unique - we've tested before at any rate??
   
   $self->invalidate_height();
   
   return scalar keys %{$self->{'_desc'}};
}


=head2 each_Descendent

 Title   : each_Descendent($sortby)
 Usage   : my @nodes = $node->each_Descendent;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : $sortby [optional] "height", "creation" or coderef to be used
           to sort the order of children nodes.

=cut

sub each_Descendent{
   my ($self, $sortby) = @_;

   # order can be based on branch length (and sub branchlength)

   $sortby ||= 'height';

   if (ref $sortby eq 'CODE') {
       return sort $sortby values %{$self->{'_desc'}};
   } else  {
       if ($sortby eq 'height') {
	   return map { $_->[0] }
		  sort { $a->[1] <=> $b->[1] || 
			 $a->[2] <=> $b->[2] } 
	       map { [$_, $_->height, $_->internal_id ] } 
	   values %{$self->{'_desc'}};
       } else {
	   return map { $_->[0] }
	          sort { $a->[1] <=> $b->[1] } 
	          map { [$_, $_->height ] }
	          values %{$self->{'_desc'}};	   
       }
   }
}

=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $node->remove_Descedent($node_foo);
 Function: Removes a specific node from being a Descendent of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects which have be previously
           passed to the add_Descendent call of this object.

=cut

sub remove_Descendent{
   my ($self,@nodes) = @_;
   foreach my $n ( @nodes ) { 
       if( $self->{'_desc'}->{$n->internal_id} ) {
	   $n->ancestor(undef);
	   $self->{'_desc'}->{$n->internal_id}->ancestor(undef);
	   delete $self->{'_desc'}->{$n->internal_id};
	   my $a1 = $self->ancestor;
	   # remove unecessary nodes if we have removed the part which branches.
	   if( $a1 ) {
	       my $bl = $self->branch_length || 0;
	       my @d = $self->each_Descendent;
	       if (scalar @d == 1) {
		   $d[0]->branch_length($bl + ($d[0]->branch_length || 0));
		   $a1->add_Descendent($d[0]);
	       }
	       $a1->remove_Descendent($self);
	   }
       
       } else { 
	   if( $self->verbose ) {
	       $self->debug(sprintf("no node %s (%s) listed as a descendent in this node %s (%s)\n",$n->id, $n,$self->id,$self));
	       $self->debug("Descendents are " . join(',', keys %{$self->{'_desc'}})."\n");
	   }
       }
   }
   1;
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
   my ($self) = @_;
   # this won't cleanup the nodes themselves if you also have
   # a copy/pointer of them (I think)...
   while( my ($node,$val) = each %{ $self->{'_desc'} } ) {
       $val->ancestor(undef);
   }
   $self->{'_desc'} = {};
   1;
}

=head2 get_all_Descendents

 Title   : get_all_Descendents
 Usage   : my @nodes = $node->get_all_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : none

=cut

# implemented in the interface 

=head2 ancestor

 Title   : ancestor
 Usage   : $obj->ancestor($newval)
 Function: Set the Ancestor
 Returns : value of ancestor
 Args    : newvalue (optional)

=cut

sub ancestor{
   my $self = shift;
   $self->{'_ancestor'} = shift @_ if @_;
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
    my $self = shift;
    $self->{'_branch_length'} = shift @_ if @_;
    return $self->{'_branch_length'};
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
    my $self = shift;
    $self->{'_bootstrap'} = shift @_ if @_;
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
    my $self = shift;
    $self->{'_description'} = shift @_ if @_;
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
    my $self = shift;
    $self->{'_id'} = shift @_ if @_;
    return $self->{'_id'};
}

sub DESTROY {
    my ($self) = @_;
    # try to insure that everything is cleaned up
    $self->SUPER::DESTROY();
    if( defined $self->{'_desc'} &&
	ref($self->{'_desc'}) =~ /ARRAY/i ) {
	while( my ($nodeid,$node) = each %{ $self->{'_desc'} } ) {
	    $node->ancestor(undef); # insure no circular references
	    $node->DESTROY();
	    $node = undef;
	}
	$self->{'_desc'} = {};
    }
}

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

sub internal_id{
   return $_[0]->_creation_id;
}


=head2 _creation_id

 Title   : _creation_id
 Usage   : $obj->_creation_id($newval)
 Function: a private method signifying the internal creation order
 Returns : value of _creation_id
 Args    : newvalue (optional)


=cut

sub _creation_id{
    my $self = shift @_;
    $self->{'_creation_id'} = shift @_ if( @_);
    return $self->{'_creation_id'} || 0;
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

sub height { 
    my ($self) = @_;

    return $self->{'_height'} if( defined $self->{'_height'} );
    
    if( $self->is_Leaf ) { 
       if( !defined $self->branch_length ) { 
	   $self->debug(sprintf("Trying to calculate height of a node when a Node (%s) has an undefined branch_length\n",$self->id || '?' ));
	   return 0;
       }
       return $self->branch_length;
   }
   my $max = 0;
   foreach my $subnode ( $self->each_Descendent ) { 
       my $s = $subnode->height;
       if( $s > $max ) { $max = $s; }
   }
   return ($self->{'_height'} = $max + ($self->branch_length || 1));
}


=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the node'e height in the tree
 Returns : nothing
 Args    : none

=cut

#'

sub invalidate_height { 
    my ($self) = @_;
    
    $self->{'_height'} = undef;
    if( $self->ancestor ) {
	$self->ancestor->invalidate_height;
    }
}

1;

