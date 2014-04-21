#
# BioPerl module for Bio::Tree::NodeI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => 'file.dnd');

    my $tree = $treeio->next_tree; # we'll assume it worked for demo purposes
                                   # you might want to test that it was defined

    my $rootnode = $tree->get_root_node;

    # process just the next generation
    foreach my $node ( $rootnode->each_Descendent() ) {
	print "branch len is ", $node->branch_length, "\n";
    }

    # process all the children
    my $example_leaf_node;
    foreach my $node ( $rootnode->get_all_Descendents() ) {
	if( $node->is_Leaf ) { 
	    print "node is a leaf ... "; 
            # for example use below
            $example_leaf_node = $node unless defined $example_leaf_node; 
	}
	print "branch len is ", $node->branch_length, "\n";
    }

    # The ancestor() method points to the parent of a node
    # A node can only have one parent

    my $parent = $example_leaf_node->ancestor;

    # parent won't likely have an description because it is an internal node
    # but child will because it is a leaf

    print "Parent id: ", $parent->id," child id: ", 
          $example_leaf_node->id, "\n";


=head1 DESCRIPTION

A NodeI is capable of the basic structure of building a tree and
storing the branch length between nodes.  The branch length is the
length of the branch between the node and its ancestor, thus a root
node in a Tree will not typically have a valid branch length.

Various implementations of NodeI may extend the basic functions and
allow storing of other information (like attatching a species object
or full sequences used to build a tree or alternative sequences).  If
you don't know how to extend a Bioperl object please ask, happy to
help, we would also greatly appreciate contributions with improvements
or extensions of the objects back to the Bioperl code base so that
others don't have to reinvent your ideas.


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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Aaron Mackey amackey@virginia.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tree::NodeI;
use strict;
no warnings 'recursion';

use base qw(Bio::Root::RootI);

=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendent($node);
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

=head2 get_all_Descendents

 Title   : get_all_Descendents($sortby)
 Usage   : my @nodes = $node->get_all_Descendents;
 Function: Recursively fetch all the nodes and their descendents
           *NOTE* This is different from each_Descendent
 Returns : Array or Bio::Tree::NodeI objects
 Args    : $sortby [optional] "height", "creation", "alpha", "revalpha", 
           or a coderef to be used to sort the order of children nodes.

=cut

sub get_all_Descendents{
   my ($self, $sortby) = @_;
   $sortby ||= 'none';   
   my @nodes;
   foreach my $node ( $self->each_Descendent($sortby) ) {
       push @nodes, ($node,$node->get_all_Descendents($sortby));
   }
   return @nodes;
}

*get_Descendents = \&get_all_Descendents;

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
       $count += 1;
       $node->can('descendent_count') ? $count += $node->descendent_count : next;
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
   return join('',defined $self->id_output ? $self->id_output : '',
		  defined $self->branch_length ? ':' . $self->branch_length 
		  : ' ')
}

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength to get to the tip.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=cut

sub height{
    my ($self) = @_;

    return 0 if( $self->is_Leaf );
    
    my $max = 0;
    foreach my $subnode ( $self->each_Descendent ) { 
	my $s = $subnode->height + $subnode->branch_length;;
	if( $s > $max ) { $max = $s; }
    }
    return $max;
}

=head2 depth

 Title   : depth
 Usage   : my $len = $node->depth
 Function: Returns the depth of the tree starting at this
           node.  Depth is the distance from this node to the root.
 Returns : The branch length to the root.
 Args    : none

=cut

sub depth{
   my ($self) = @_;
   
   my $depth = 0;
   my $node = $self;
   while( defined $node->ancestor ) { 
       $depth += $node->branch_length;
       $node = $node->ancestor;
   }
   return $depth;
}

=head2 Get/Set methods

=cut

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: Get/Set the branch length
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
 Function: The human readable identifier for the node 
 Returns : value of human readable id
 Args    : newvalue (optional)


=cut

sub id{
    my ($self)= @_;
    $self->throw_not_implemented();
}

=head2 internal_id

 Title   : internal_id
 Usage   : my $internalid = $node->internal_id
 Function: Returns the internal unique id for this Node
 Returns : unique id
 Args    : none

=cut

sub internal_id{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: Get/Set the description string
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
 Function: Get/Set the bootstrap value
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
 Function: Get/Set the ancestor node pointer for a Node
 Returns : Null if this is top level node
 Args    : none

=cut


sub ancestor{
   my ($self,@args) = @_;
    $self->throw_not_implemented();
}

=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the node height in the tree
 Returns : nothing
 Args    : none

=cut

sub invalidate_height { 
    shift->throw_not_implemented();
}

=head2 Methods for associating Tag/Values with a Node

These methods associate tag/value pairs with a Node

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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $node->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag {
    shift->throw_not_implemented();
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $node->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None


=cut

sub remove_all_tags{
    shift->throw_not_implemented();  
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $node->get_all_tags()
 Function: Gets all the tag names for this Node
 Returns : Array of tagnames
 Args    : None


=cut

sub get_all_tags {
    shift->throw_not_implemented();
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $node->get_tag_values($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name


=cut

sub get_tag_values{
    shift->throw_not_implemented();
}

=head2 has_tag

 Title   : has_tag
 Usage   : $node->has_tag($tag)
 Function: Boolean test if tag exists in the Node
 Returns : Boolean
 Args    : $tag - tagname


=cut

sub has_tag{
    shift->throw_not_implemented();
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

1;
