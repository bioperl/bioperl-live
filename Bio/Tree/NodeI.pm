# NodeI.pm,v 1.10 2002/04/18 12:52:46 jason Exp
#
# BioPerl module for Bio::Tree::NodeI
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::NodeI - Interface describing a Tree Node

=head1 SYNOPSIS

    # get a Tree::NodeI somehow
    # like from a TreeIO
    use Bio::TreeIO;
    # read in a clustalw NJ in phylip/newick format
    my $treeio = new Bio::TreeIO(-format => 'newick', -file => 'file.dnd');
    my $tree = $treeio->next_tree; # we'll assume it worked for demo purposes

    my $rootnode = $tree->get_root_node;

    # process just the next generation
    foreach my $node ( $tree->each_Descendent() ) {
	print "branch len is ", $node->branch_length, "\n";
    }
    # process all the children
    foreach my $node ( $tree->get_Descendents() ) {
	if( $node->is_Leaf ) { 
	    print "node is a leaf ... "; 
	}
	print "branch len is ", $node->branch_length, "\n";
    }    

=head1 DESCRIPTION

A NodeI is capable of the basic structure of building a tree and
storing the branch length between nodes.  The branch length is the
length of the branch between the node and its ancestor, thus a root
node in a Tree will not typically have a valid branch length.

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

package Bio::Tree::NodeI;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;
@ISA = qw(Bio::Root::RootI);


=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendant($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI


=cut

sub add_Descendent{
   my ($self,@args) = @_;

   $self->throw_not_implemented();
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
   $self->throw_not_implemented();
}

=head2 Decorated Interface methods

=cut

=head2 get_Descendents

 Title   : get_Descendents
 Usage   : my @nodes = $node->get_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : none

=cut

sub get_Descendents{
   my ($self) = @_;
   my @nodes;
   foreach my $node ( $self->each_Descendent ) {
       push @nodes, ($node->get_Descendents, $node);
   }
   return @nodes;
}

=head2 is_Leaf

 Title   : is_Leaf
 Usage   : if( $node->is_Leaf ) 
 Function: Get Leaf status
 Returns : boolean
 Args    : none

=cut

sub is_Leaf{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 descendent_count

 Title   : descendent_count
 Usage   : my $count = $node->descendent_count;
 Function: Counts the number of descendents a node has 
           (and all of their subnodes)
 Returns : integer
 Args    : none

=cut

sub descendent_count{
   my ($self) = @_;
   my $count = 0;
   
   foreach my $node ( $self->each_Descendent ) { 
       $count += 1 + $node->descendent_count;
   }
   return $count;
}

=head2 to_string

 Title   : to_string
 Usage   : my $str = $node->to_string()
 Function: For debugging, provide a node as a string
 Returns : string
 Args    : none


=cut

sub to_string{
   my ($self) = @_;
   return sprintf("%s%s",defined $self->id ? $self->id : '',
		  defined $self->branch_length ? 
		  ":" . $self->branch_length : ' ');
}

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

sub height{
   my ($self) = @_;
   
   if( $self->is_Leaf ) { 
       if( !defined $self->branch_length ) { 
	   $self->debug(sprintf("Trying to calculate height of a node when a Node (%s) has an undefined branch_length",$self->id || '?' ));
	   return 0;
       }
       return $self->branch_length;
   }
   my $max = 0;
   foreach my $subnode ( $self->each_Descendent ) { 
       my $s = $subnode->height;
       if( $s > $max ) { $max = $s; }
   }
   return $max + ($self->branch_length || 0);
}

=head2 Get/Set methods

=cut

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: 
 Example : 
 Returns : value of branch_length
 Args    : newvalue (optional)


=cut

sub branch_length{
    my ($self)= @_;
    $self->throw_not_implemented();
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
    my ($self)= @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
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
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 ancestor

 Title   : ancestor
 Usage   : my $node = $node->ancestor;
 Function: Get/Set a Node's ancestor node
 Returns : Null if this is top level node
 Args    : none

=cut

#'
sub ancestor{
   my ($self,@args) = @_;
    $self->throw_not_implemented();
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
    shift->throw_not_implemented();
}

1;
