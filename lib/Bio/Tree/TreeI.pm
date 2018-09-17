#
# BioPerl module for Bio::Tree::TreeI
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

Bio::Tree::TreeI - A Tree object suitable for lots of things, designed
  originally for Phylogenetic Trees.

=head1 SYNOPSIS

  # get a Bio::Tree::TreeI somehow
  # like from a TreeIO
  my $treeio = Bio::TreeIO->new(-format => 'newick', -file => 'treefile.dnd');
  my $tree   = $treeio->next_tree;
  my @nodes  = $tree->get_nodes;
  my @leaves = $tree->get_leaf_nodes;
  my $root   = $tree->get_root_node;

=head1 DESCRIPTION

This object holds a pointer to the Root of a Tree which is a
Bio::Tree::NodeI.

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

Aaron Mackey, amackey@virginia.edu
Elia Stupka,  elia@fugu-sg.org
Sendu Bala,   bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::TreeI;
use strict;

use base qw(Bio::Root::RootI);

=head2 get_nodes

 Title   : get_nodes
 Usage   : my @nodes = $tree->get_nodes()
 Function: Return list of Tree::NodeI objects
 Returns : array of Tree::NodeI objects
 Args    : (named values) hash with one value 
           order => 'b|breadth' first order or 'd|depth' first order

=cut

sub get_nodes{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_root_node

 Title   : get_root_node
 Usage   : my $node = $tree->get_root_node();
 Function: Get the Top Node in the tree, in this implementation
           Trees only have one top node.
 Returns : Bio::Tree::NodeI object
 Args    : none

=cut

sub get_root_node{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 number_nodes

 Title   : number_nodes
 Usage   : my $size = $tree->number_nodes
 Function: Find the number of nodes in the tree.
 Returns : int
 Args    : none

=cut

sub number_nodes{
   my ($self) = @_;
   my $root = $self->get_root_node;
   if( defined $root && $root->isa('Bio::Tree::NodeI'))  {
       return ($root->descendent_count + 1);
   }
   return 0;
}

=head2 total_branch_length

 Title   : total_branch_length
 Usage   : my $size = $tree->total_branch_length
 Function: Returns the sum of the length of all branches
 Returns : integer
 Args    : none

=cut

sub total_branch_length {
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 height

 Title   : height
 Usage   : my $height = $tree->height
 Function: Gets the height of tree - this LOG_2($number_nodes)
           WARNING: this is only true for strict binary trees.  The TreeIO
           system is capable of building non-binary trees, for which this
           method will currently return an incorrect value!!
 Returns : integer
 Args    : none

=cut

sub height{
   my ($self) = @_;
   my $nodect =  $self->number_nodes;
   return 0 if( ! $nodect ); 
   return log($nodect) / log(2);
}

=head2 id

 Title   : id
 Usage   : my $id = $tree->id();
 Function: An id value for the tree
 Returns : scalar
 Args    : 


=cut

sub id{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 score

 Title   : score
 Usage   : $obj->score($newval)
 Function: Sets the associated score with this tree
           This is a generic slot which is probably best used 
           for log likelihood or other overall tree score
 Returns : value of score
 Args    : newvalue (optional)


=cut

sub score{
   my ($self,$value) = @_;
   $self->throw_not_implemented();
}

=head2 get_leaf_nodes

 Title   : get_leaf_nodes
 Usage   : my @leaves = $tree->get_leaf_nodes()
 Function: Returns the leaves (tips) of the tree
 Returns : Array of Bio::Tree::NodeI objects
 Args    : none


=cut

sub get_leaf_nodes{
   my ($self) = @_;
   return grep { $_->is_Leaf() } $self->get_nodes(-sortby  => 'none');
}


=head2 Methods for associating Tag/Values with a Tree

These methods associate tag/value pairs with a Tree

=head2 set_tag_value

 Title   : set_tag_value
 Usage   : $tree->set_tag_value($tag,$value)
           $tree->set_tag_value($tag,@values)
 Function: Sets a tag value(s) to a tree. Replaces old values.
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub set_tag_value{
    shift->throw_not_implemented();
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $tree->add_tag_value($tag,$value)
 Function: Adds a tag value to a tree 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag


=cut

sub add_tag_value{
    shift->throw_not_implemented();
}

=head2 remove_tag

 Title   : remove_tag
 Usage   : $tree->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove


=cut

sub remove_tag {
    shift->throw_not_implemented();
}

=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $tree->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None


=cut

sub remove_all_tags{
    shift->throw_not_implemented();
}

=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $tree->get_all_tags()
 Function: Gets all the tag names for this Tree
 Returns : Array of tagnames
 Args    : None


=cut

sub get_all_tags {
    shift->throw_not_implemented();
}

=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $tree->get_tag_values($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name


=cut

sub get_tag_values{
    shift->throw_not_implemented();
}

=head2 has_tag

 Title   : has_tag
 Usage   : $tree->has_tag($tag)
 Function: Boolean test if tag exists in the Tree
 Returns : Boolean
 Args    : $tag - tagname


=cut

sub has_tag{
    shift->throw_not_implemented();
}

1;
