#
# BioPerl module for Bio::Tree::Tree
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

Bio::Tree::Tree - An implementation of the TreeI interface.

=head1 SYNOPSIS

    use Bio::TreeIO;

    # like from a TreeIO
    my $treeio = Bio::TreeIO->new(-format => 'newick', -file => 'treefile.dnd');
    my $tree = $treeio->next_tree;
    my @nodes = $tree->get_nodes;
    my $root = $tree->get_root_node;

=head1 DESCRIPTION

This object holds handles to Nodes which make up a tree.

=head1 IMPLEMENTATION NOTE

This implementation of Bio::Tree::Tree contains Bio::Tree:::NodeI; mainly linked
via the root node. As NodeI can potentially contain circular references (as
nodes will need to refer to both parent and child nodes), Bio::Tree::Tree will
remove those circular references when the object is garbage-collected. This has
some side effects; primarily, one must keep the Tree in scope or have at least
one reference to it if working with nodes. The fix is to count the references to
the nodes and if it is greater than expected retain all of them, but it requires
an additional prereq and thus may not be worth the effort.  This only shows up
in minor edge cases, though (see Bug #2869).

Example of issue:

  # tree is not assigned to a variable, so passes from memory after
  # root node is passed
  my $root = Bio::TreeIO->new(-format => 'newick', -file => 'foo.txt')->next_tree
                 ->get_root_node;

  # gets nothing, as all Node links are broken when Tree is garbage-collected above
  my @descendents = $root->get_all_Descendents;

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
Sendu Bala   bix@sendu.me.uk
Mark A. Jensen maj@fortinbras.us

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::Tree;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::Tree::TreeI Bio::Tree::TreeFunctionsI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Tree->new();
 Function: Builds a new Bio::Tree::Tree object 
 Returns : Bio::Tree::Tree
 Args    : -root     => L<Bio::Tree::NodeI> object which is the root
              OR
           -node     => L<Bio::Tree::NodeI> object from which the root will be
                        determined
 
           -nodelete => boolean, whether or not to try and cleanup all
                        the nodes when this this tree goes out of scope.
           -id       => optional tree ID
           -score    => optional tree score value

=cut

sub new {
    my ($class, @args) = @_;
  
    my $self = $class->SUPER::new(@args);
    $self->{'_rootnode'} = undef;
    $self->{'_maxbranchlen'} = 0;
    $self->_register_for_cleanup(\&cleanup_tree);
    my ($root, $node, $nodel, $id, $score) =
        $self->_rearrange([qw(ROOT NODE NODELETE ID SCORE)], @args);
  
    if ($node && ! $root) {
        $self->throw("Must supply a Bio::Tree::NodeI") unless ref($node) && $node->isa('Bio::Tree::NodeI');
        my @lineage = $self->get_lineage_nodes($node);
        $root = shift(@lineage) || $node;
    
        # to stop us pulling in entire database of a Bio::Taxon when we later do
        # get_nodes() or similar, specifically set ancestor() for each node
        if ($node->isa('Bio::Taxon')) {
            push(@lineage, $node) unless $node eq $root;
            my $ancestor = $root;
            foreach my $lineage_node (@lineage) {
                $lineage_node->ancestor($ancestor);
            } continue { $ancestor = $lineage_node; }
        }
    }
    if ($root) {
        $self->set_root_node($root);
    }
  
    $self->nodelete($nodel || 0);
    $self->id($id)       if defined $id;
    $self->score($score) if defined $score;
    return $self;
}


=head2 nodelete

 Title   : nodelete
 Usage   : $obj->nodelete($newval)
 Function: Get/Set Boolean whether or not to delete the underlying
           nodes when it goes out of scope.  By default this is false
           meaning trees are cleaned up.
 Returns : boolean
 Args    : on set, new boolean value

=cut

sub nodelete {
    my $self = shift;
    return $self->{'nodelete'} = shift if @_;
    return $self->{'nodelete'};
}


=head2 get_nodes

 Title   : get_nodes
 Usage   : my @nodes = $tree->get_nodes()
 Function: Return list of Bio::Tree::NodeI objects
 Returns : array of Bio::Tree::NodeI objects
 Args    : (named values) hash with one value 
           order => 'b|breadth' first order or 'd|depth' first order
           sortby => [optional] "height", "creation", "alpha", "revalpha",
           or coderef to be used to sort the order of children nodes. See L<Bio::Tree::Node> for details

=cut

sub get_nodes {
    my ($self, @args) = @_;
    my ($order, $sortby) = $self->_rearrange([qw(ORDER SORTBY)], @args);
    $order  ||= 'depth';
    $sortby ||= 'none';

    my @children;
    my $node = $self->get_root_node;
    if ($node) {
        if ($order =~ m/^b/oi) { # breadth-first
            @children = ($node);
            my @to_process = ($node);
            while( @to_process ) { 
                my $n = shift @to_process;
                my @c  = $n->each_Descendent($sortby);
                push @children, @c;
                push @to_process, @c;
            }
        } elsif ($order =~ m/^d/oi) { # depth-first
            @children = ($node, $node->get_all_Descendents($sortby));
        } else {
            $self->verbose(1);
            $self->warn("specified an order '$order' which I don't understan\n");
        }
    }

    return @children;
}


=head2 get_root_node

 Title   : get_root_node
 Usage   : my $node = $tree->get_root_node();
 Function: Get the Top Node in the tree, in this implementation
           Trees only have one top node.
 Returns : Bio::Tree::NodeI object
 Args    : none

=cut

sub get_root_node {
    my ($self) = @_;
    return $self->{'_rootnode'};
}


=head2 set_root_node

 Title   : set_root_node
 Usage   : $tree->set_root_node($node)
 Function: Set the Root Node for the Tree
 Returns : Bio::Tree::NodeI
 Args    : Bio::Tree::NodeI

=cut

sub set_root_node {
    my $self = shift;
    if ( @_ ) { 
        my $value = shift;
        if ( defined $value && ! $value->isa('Bio::Tree::NodeI') ) { 
            $self->warn("Trying to set the root node to $value which is not a Bio::Tree::NodeI");
            return $self->get_root_node;
        }
        $self->{'_rootnode'} = $value;
    } 
    return $self->get_root_node;
}


=head2 total_branch_length

 Title   : total_branch_length
 Usage   : my $size = $tree->total_branch_length
 Function: Returns the sum of the length of all branches
 Returns : real
 Args    : none

=cut

sub total_branch_length { shift->subtree_length }


=head2 subtree_length

 Title   : subtree_length
 Usage   : my $subtree_size = $tree->subtree_length($internal_node)
 Function: Returns the sum of the length of all branches in a subtree
           under the node. Calculates the size of the whole tree
           without an argument (but only if root node is defined)
 Returns : real or undef
 Args    : Bio::Tree::NodeI object, defaults to the root node

=cut

sub subtree_length {
    my $tree = shift;
    my $node = shift || $tree->get_root_node;
    return unless $node;
    my $sum = 0;
    for ( $node->get_all_Descendents ) {
        $sum += $_->branch_length || 0;
    }
    return $sum;
}


=head2 id

 Title   : id
 Usage   : my $id = $tree->id();
 Function: An id value for the tree
 Returns : scalar
 Args    : [optional] new value to set

=cut

sub id {
   my ($self, $val) = @_;
   if ( defined $val ) { 
       $self->{'_treeid'} = $val;
   }
   return $self->{'_treeid'};
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

sub score {
   my ($self, $val) = @_;
   if ( defined $val ) { 
       $self->{'_score'} = $val;
   }
   return $self->{'_score'};
}


# decorated interface TreeI Implements this

=head2 height

 Title   : height
 Usage   : my $height = $tree->height
 Function: Gets the height of tree - this LOG_2($number_nodes)
           WARNING: this is only true for strict binary trees.  The TreeIO
           system is capable of building non-binary trees, for which this
           method will currently return an incorrect value!!
 Returns : integer
 Args    : none

=head2 number_nodes

 Title   : number_nodes
 Usage   : my $size = $tree->number_nodes
 Function: Returns the number of nodes in the tree
 Returns : integer
 Args    : none

=head2 as_text

 Title   : as_text
 Usage   : my $tree_as_string = $tree->as_text($format)
 Function: Returns the tree as a string representation in the 
           desired format, e.g.: 'newick', 'nhx' or 'tabtree' (the default)
 Returns : scalar string
 Args    : format type as specified by Bio::TreeIO
 Note    : This method loads the Bio::TreeIO::$format module
           on the fly, and commandeers the _write_tree_Helper
           routine therein to create the tree string. 

=cut

sub as_text {
    my $self = shift;
    my $format = shift || 'tabtree';
    my $params_input = shift || {};

    my $iomod = "Bio::TreeIO::$format";
    $self->_load_module($iomod);

    my $string = '';
    open my $fh, '>', \$string or $self->throw("Could not write '$string' as file: $!");
    my $test = $iomod->new( -format => $format, -fh => $fh );

    # Get the default params for the given IO module.
    $test->set_params($params_input);

    $test->write_tree($self);
    close $fh;
    return $string;
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

sub set_tag_value {
    my ($self, $tag, @values) = @_;
    if ( ! defined $tag || ! scalar @values ) {
        $self->warn("cannot call set_tag_value with an undefined value");
    }
    $self->remove_tag ($tag);
    map { push @{$self->{'_tags'}->{$tag}}, $_ } @values;
    return scalar @{$self->{'_tags'}->{$tag}};
}


=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $tree->add_tag_value($tag,$value)
 Function: Adds a tag value to a tree 
 Returns : number of values stored for this tag
 Args    : $tag   - tag name
           $value - value to store for the tag

=cut

sub add_tag_value {
    my ($self, $tag, $value) = @_;
    if ( ! defined $tag || ! defined $value ) {
        $self->warn("cannot call add_tag_value with an undefined value");
    }
    push @{$self->{'_tags'}->{$tag}}, $value;
    return scalar @{$self->{'_tags'}->{$tag}};
}


=head2 remove_tag

 Title   : remove_tag
 Usage   : $tree->remove_tag($tag)
 Function: Remove the tag and all values for this tag
 Returns : boolean representing success (0 if tag does not exist)
 Args    : $tag - tagname to remove

=cut

sub remove_tag {
    my ($self, $tag) = @_;
    if ( exists $self->{'_tags'}->{$tag} ) {
        $self->{'_tags'}->{$tag} = undef;
        delete $self->{'_tags'}->{$tag};
        return 1;
    }
    return 0;
}


=head2 remove_all_tags

 Title   : remove_all_tags
 Usage   : $tree->remove_all_tags()
 Function: Removes all tags 
 Returns : None
 Args    : None

=cut

sub remove_all_tags {
    my ($self) = @_;
    $self->{'_tags'} = {};
    return;
}


=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $tree->get_all_tags()
 Function: Gets all the tag names for this Tree
 Returns : Array of tagnames
 Args    : None

=cut

sub get_all_tags {
    my ($self) = @_;
    my @tags = sort keys %{$self->{'_tags'} || {}};
    return @tags;
}


=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $tree->get_tag_values($tag)
 Function: Gets the values for given tag ($tag)
 Returns : Array of values or empty list if tag does not exist
 Args    : $tag - tag name

=cut

sub get_tag_values {
    my ($self, $tag) = @_;
    return wantarray ? @{$self->{'_tags'}->{$tag} || []} :
                      (@{$self->{'_tags'}->{$tag} || []})[0];
}


=head2 has_tag

 Title   : has_tag
 Usage   : $tree->has_tag($tag)
 Function: Boolean test if tag exists in the Tree
 Returns : Boolean
 Args    : $tag - tagname

=cut

sub has_tag {
    my ($self, $tag) = @_;
    return exists $self->{'_tags'}->{$tag};
}


# safe tree clone that doesn't seg fault

=head2 clone

 Title   : clone
 Alias   : _clone
 Usage   : $tree_copy = $tree->clone();
           $subtree_copy = $tree->clone($internal_node);
 Function: Safe tree clone that doesn't segfault
 Returns : Bio::Tree::Tree object
 Args    : [optional] $start_node, Bio::Tree::Node object

=cut

sub clone {
    my ($self, $parent, $parent_clone) = @_;
    $parent ||= $self->get_root_node;
    $parent_clone ||= $self->_clone_node($parent);

    foreach my $node ($parent->each_Descendent()) {
        my $child = $self->_clone_node($node);
        $child->ancestor($parent_clone);
        $self->_clone($node, $child);
    }
    $parent->ancestor && return;

    my $tree = $self->new(-root => $parent_clone);
    return $tree;
}


# -- private internal methods --

sub cleanup_tree {
    my $self = shift;
    unless( $self->nodelete ) {
        for my $node ($self->get_nodes(-order  => 'b', -sortby => 'none')) {
            $node->node_cleanup;
        }
    }
    $self->{'_rootnode'} = undef;
}

1;
