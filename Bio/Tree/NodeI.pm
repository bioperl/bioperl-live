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

=head1 DESCRIPTION

The NodeI interface describes the core methods necessary to model a
basic rooted tree structure as a set of Perl objects.

An implementing class need only implement the five methods in NodeI
(parent, children, add_child, remove_child, branch_length) to gain all
of the functionality contained within NodeFunctionsI, which is a set
of object methods built to perform more advanced tree manipulations
using only the basic methods of the NodeI interface. An implementing
class should use both Bio::Tree::NodeI and Bio::Tree::NodeFunctionsI
as a base.

Note that the only method in NodeI which is a combined getter / setter
method is branch_length; all the other methods either query the tree
structure (parent / children) to return an object or array, or perform
an action on the tree structure (add_child / remove_child) where the
return value doesn't matter.

Nearly all of the code for performing tree manipulations is contained
within NodeFunctionsI. Although there are similarly named TreeI,
TreeFunctionsI and Tree classes to parallel the Node classes, a Tree
object simply contains a pointer to the root node and some additional
metadata, such as score, id, and description. Most Tree methods that
perform tree manipulations simply pass the call on to the root
node. See the TreeI, TreeFunctionsI and Tree modules for more
information.

Various implementations of NodeI may extend the basic functions and
allow storing of other information (like attatching a species object
or full sequences used to build a tree or alternative sequences). For
example, most phylogenetic trees store labels and bootstrap values; in
the stock BioPerl implementation (Node.pm) the method node.id() is
implemented to fulfill the LabeledNodeI.pm interface, while the
node.bootstrap() method is unique to the Node.pm implementation and
relies on the get_tag_value and set_tag_value methods provided by the
Bio::Tree::TagValueHolder base class.

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

  https://redmine.open-bio.org/projects/bioperl/

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

=head2 parent

 Title   : parent
 Usage   : my $parent = $node->parent;
 Function: Return this node's parent, or undef if this node has no parent.
 Returns : Bio::Tree::NodeI or undef if no parent exists (i.e. the
           current node is the root)
 Args    : none

=cut

sub parent{
    shift->throw_not_implemented();
}
sub ancestor { shift->parent(@_) }

=head2 branch_length

 Title   : branch_length
 Usage   : $node->branch_length(0.5);
 Function: Get/Set the length of the branch between the current node and its parent
 Returns : Branch length
 Args    : new branch length (optional)

=cut

sub branch_length{
    shift->throw_not_implemented();
}
sub distance_to_parent { shift->branch_length(@_) }

=head2 children

 Title   : children
 Usage   : my @nodes = $node->children;

 Function: List the direct children of this node (but not all
           descendents; use .nodes() or .leaves() for this). The order
           of children need not be defined or stable in order to
           implement this interface, but the stock Bioperl
           implementation (in L<Bio::Tree::Node::children>) sorts
           children by an internally-stored sort value.

 Returns : Array of Bio::Tree::NodeI objects
 Args    : none

=cut

sub children{
   shift->throw_not_implemented();
}
sub each_Descendent { shift->children(@_)}

=head2 add_child

 Title   : add_child
 Usage   : $node->add_child($new_node);
 Function: Add a child to a node.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the node to add as a child of the current node

=cut

sub add_child{
   shift->throw_not_implemented();
}
sub add_Descendent { shift->add_child(@_) }

=head2 remove_child

 Title   : remove_child
 Usage   : $node->remove_child($child_node);
 Function: Removes a child from the current node.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the child to remove from the current node

=cut

sub remove_child {
    shift->throw_not_implemented();
}
sub remove_Descendent { shift->remove_child(@_) }

=head2 remove_child

 Title   : remove_child
 Usage   : $node->remove_child($child_node);
 Function: Removes a child from the current node.
 Returns : Nothing
 Args    : Bio::Tree::NodeI, the child to remove from the current node

=cut

sub tree {
    shift->throw_not_implemented();
}

1;
