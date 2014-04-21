#
# BioPerl module for Bio::Tree::AlleleNode
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

Bio::Tree::AlleleNode - A Node with Alleles attached

=head1 SYNOPSIS

  use Bio::Tree::AlleleNode;

=head1 DESCRIPTION

AlleleNodes are basic L<Bio::Tree::Node>s with the added ability to
add Genotypes alleles as defined by the L<Bio::PopGen::IndividualI>
interface.  Genotypes are defined by the L<Bio::PopGen::GenotypeI>
interface, you will probably want to use the L<Bio::PopGen::Genotype>
implementation.

This is implemented via containment to avoid multiple inheritance
problems.  Their is a L<Bio::PopGen::Individual> object which handles
the L<Bio::PopGen::IndividualI> interface, and is accessible via the
L<Bio::Tree::AlleleNode::individual> method.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=head1 HISTORY

This module was re-written to be a combination of
L<Bio::PopGen::Individual> and L<Bio::Tree::Node> primarily for use in
L<Bio::PopGen::Simulation::Coalescent> simulations.

=cut

# Let the code begin...


package Bio::Tree::AlleleNode;
use vars qw($UIDCOUNTER);
use strict;
BEGIN { $UIDCOUNTER = 1 }

use Bio::PopGen::Individual;
use Bio::PopGen::Genotype;

use base qw(Bio::Tree::Node Bio::PopGen::IndividualI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::AlleleNode->new();
 Function: Builds a new Bio::Tree::AlleleNode() object 
 Returns : an instance of Bio::Tree::AlleleNode
 Args    : -unique_id     => $id,
           -genotypes     => \@genotypes
           -left          => pointer to Left descendent (optional)
           -right         => pointer to Right descenent (optional)
	   -branch_length => branch length [integer] (optional)
           -bootstrap     => value   bootstrap value (string)
           -description   => description of node
           -id            => human readable (unique) id for node
                             Should NOT contain the characters 
                             '();:'
=cut

sub new { 
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    $self->individual( Bio::PopGen::Individual->new(@args));
    return $self;
}

=head2 individual

 Title   : individual
 Usage   : $obj->individual($newval)
 Function: Get/Set Access to the underlying individual object
 Returns : L<Bio::PopGen::Individual> object
 Args    : on set, new value (L<Bio::PopGen::Individual>)


=cut

sub individual {
    my ($self,$newval) = @_;
    if( defined $newval || ! defined $self->{'individual'} ) {
	$newval = Bio::PopGen::Individual->new() unless defined $newval;
	$self->{'individual'} = $newval;
    }
    return $self->{'individual'};
}

=head2 Bio::PopGen::Individual methods

Methods required by L<Bio::PopGen::IndividualI>.


=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $individual->unique_id
 Function: Unique Identifier
 Returns : string representing unique identifier
 Args    : string


=cut

sub unique_id{
    my $self = shift;
    $self->individual->unique_id(@_);
}

=head2 num_of_results

 Title   : num_of_results
 Usage   : my $count = $person->num_results;
 Function: returns the count of the number of Results for a person
 Returns : integer
 Args    : none

=cut

sub num_of_results {
    my $self = shift;
    $self->individual->num_of_results(@_);
}

=head2 add_Genotype

 Title   : add_Genotype
 Usage   : $individual->add_Genotype
 Function: add a genotype value, only a single genotype
           may be associated 
 Returns : count of the number of genotypes associated with this individual
 Args    : @genotypes - Bio::PopGen::GenotypeI object(s) containing 
                        alleles plus a marker name

=cut

sub add_Genotype {
    my $self = shift;
    $self->individual->add_Genotype(@_);
}

=head2 reset_Genotypes

 Title   : reset_Genotypes
 Usage   : $individual->reset_Genotypes;
 Function: Reset the genotypes stored for this individual
 Returns : none
 Args    : none


=cut

sub reset_Genotypes{
    my $self = shift;
    $self->individual->reset_Genotypes(@_);
}

=head2 remove_Genotype

 Title   : remove_Genotype
 Usage   : $individual->remove_Genotype(@names)
 Function: Removes the genotypes for the requested markers
 Returns : none
 Args    : Names of markers 


=cut

sub remove_Genotype{
    my $self = shift;
    $self->individual->remove_Genotype(@_);
}

=head2 get_Genotypes

 Title   : get_Genotypes
 Usage   : my @genotypes = $ind->get_Genotypes(-marker => $markername);
 Function: Get the genotypes for an individual, based on a criteria
 Returns : Array of genotypes
 Args    : either none (return all genotypes) or 
           -marker => name of marker to return (exact match, case matters)


=cut

sub get_Genotypes{
    my $self = shift;
    $self->individual->get_Genotypes(@_);
}

=head2 has_Marker

 Title   : has_Marker
 Usage   : if( $ind->has_Marker($name) ) {}
 Function: Boolean test to see if an Individual has a genotype 
           for a specific marker
 Returns : Boolean (true or false)
 Args    : String representing a marker name


=cut

sub has_Marker{
    my $self = shift;
    $self->individual->has_Marker(@_);
}

=head2 get_marker_names

 Title   : get_marker_names
 Usage   : my @names = $individual->get_marker_names;
 Function: Returns the list of known marker names
 Returns : List of strings
 Args    : none


=cut

sub get_marker_names{
   my $self = shift;
   $self->individual->get_marker_names(@_);
}

=head2 Bio::Tree::Node methods

Methods inherited from L<Bio::Tree::Node>.


=head2 add_Descendent

 Title   : add_Descendent
 Usage   : $node->add_Descendent($node);
 Function: Adds a descendent to a node
 Returns : number of current descendents for this node
 Args    : Bio::Node::NodeI
           boolean flag, true if you want to ignore the fact that you are
           adding a second node with the same unique id (typically memory 
           location reference in this implementation).  default is false and 
           will throw an error if you try and overwrite an existing node.


=head2 each_Descendent

 Title   : each_Descendent($sortby)
 Usage   : my @nodes = $node->each_Descendent;
 Function: all the descendents for this Node (but not their descendents
					      i.e. not a recursive fetchall)
 Returns : Array of Bio::Tree::NodeI objects
 Args    : $sortby [optional] "height", "creation" or coderef to be used
           to sort the order of children nodes.


=head2 remove_Descendent

 Title   : remove_Descendent
 Usage   : $node->remove_Descedent($node_foo);
 Function: Removes a specific node from being a Descendent of this node
 Returns : nothing
 Args    : An array of Bio::Node::NodeI objects which have be previously
           passed to the add_Descendent call of this object.


=head2 remove_all_Descendents

 Title   : remove_all_Descendents
 Usage   : $node->remove_All_Descendents()
 Function: Cleanup the node's reference to descendents and reset
           their ancestor pointers to undef, if you don't have a reference
           to these objects after this call they will be cleaned up - so
           a get_nodes from the Tree object would be a safe thing to do first
 Returns : nothing
 Args    : none



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


=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length()
 Function: Get/Set the branch length
 Returns : value of branch_length
 Args    : newvalue (optional)


=head2 bootstrap

 Title   : bootstrap
 Usage   : $obj->bootstrap($newval)
 Function: Get/Set the bootstrap value
 Returns : value of bootstrap
 Args    : newvalue (optional)


=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: Get/Set the description string
 Returns : value of description
 Args    : newvalue (optional)


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: The human readable identifier for the node 
 Returns : value of human readable id
 Args    : newvalue (optional)
 Note    : id cannot contain the chracters '();:'

"A name can be any string of printable characters except blanks,
colons, semicolons, parentheses, and square brackets. Because you may
want to include a blank in a name, it is assumed that an underscore
character ("_") stands for a blank; any of these in a name will be
converted to a blank when it is read in."

from L<http://evolution.genetics.washington.edu/phylip/newicktree.html>

=cut

=head2 internal_id

 Title   : internal_id
 Usage   : my $internalid = $node->internal_id
 Function: Returns the internal unique id for this Node
           (a monotonically increasing number for this in-memory implementation
            but could be a database determined unique id in other 
	    implementations)
 Returns : unique id
 Args    : none


=head2 Bio::Node::NodeI decorated interface implemented

The following methods are implemented by L<Bio::Node::NodeI> decorated
interface.

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

=head2 height

 Title   : height
 Usage   : my $len = $node->height
 Function: Returns the height of the tree starting at this
           node.  Height is the maximum branchlength.
 Returns : The longest length (weighting branches with branch_length) to a leaf
 Args    : none

=head2 invalidate_height

 Title   : invalidate_height
 Usage   : private helper method
 Function: Invalidate our cached value of the node's height in the tree
 Returns : nothing
 Args    : none

=cut

#'

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $node->add_tag_value($tag,$value)
 Function: Adds a tag value to a node 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag


=head2 remove_tag

 Title   : remove_tag
 Usage   : $node->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove



=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $node->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None



=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $node->get_all_tags()
 Function: Gets all the tag names for this Node
 Returns : Array of tagnames
 Args    : None


=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $node->get_tag_value($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name


=head2 has_tag

 Title   : has_tag
 Usage   : $node->has_tag($tag)
 Function: Boolean test if tag exists in the Node
 Returns : Boolean
 Args    : $tag - tagname


=cut


1;
