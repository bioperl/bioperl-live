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

  http://bugzilla.open-bio.org/

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


sub parent{
    shift->throw_not_implemented();
}
sub ancestor { shift->parent(@_) }


sub add_child{
   shift->throw_not_implemented();
}
sub add_Descendent { shift->add_child(@_) }


sub remove_child {
    shift->throw_not_implemented();
}
sub remove_Descendent { shift->remove_child(@_) }


sub children{
   shift->throw_not_implemented();
}
sub each_Descendent { shift->children(@_)}


sub set_child_order {
    shift->throw_not_implemented();
}


sub is_leaf{
    shift->throw_not_implemented();
}
sub is_Leaf { shift->is_leaf(@_) }


sub branch_length{
    shift->throw_not_implemented();
}
sub distance_to_parent { shift->branch_length(@_) }


sub id {
    shift->throw_not_implemented();
}
sub name { shift->id(@_) }
sub label { shift->label(@_) }


sub internal_id {
   shift->throw_not_implemented();
}
sub unique_id { shift->internal_id(@_) }


sub bootstrap{
    shift->throw_not_implemented();
}

1;
