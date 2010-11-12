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

Bio::Tree::Tree - An Implementation of TreeI interface.

=head1 SYNOPSIS

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

  http://bugzilla.open-bio.org/

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

use base qw(Bio::Root::Root Bio::Tree::TreeI Bio::Tree::TreeFunctionsI Bio::Tree::TagValueHolder);

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
                                 the nodes when this this tree goes out
                                 of scope.
           -id       => optional tree ID
           -score    => optional tree score value

=cut

sub new {
  my($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args);
  $self->{'_rootnode'} = undef;
  $self->{'_maxbranchlen'} = 0;
  $self->_register_for_cleanup(\&cleanup_tree);
  my ($root,$node,$nodel,$id,$score)= $self->_rearrange([qw(ROOT NODE NODELETE 
                              ID SCORE)], @args);
  
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

sub rooted {
    my $self = shift;
    $self->{'_rooted'} = shift if @_;
    if (defined $self->{'_rooted'}) {
	return $self->{'_rooted'};
    } else {
	# Default to a rooted tree.
	return 1;	
    }
}
sub set_rooted { shift->rooted(@_) }
sub is_rooted { shift->rooted }

sub nodelete{
    my $self = shift;
    return $self->{'nodelete'} = shift if @_;
    return $self->{'nodelete'};
}

sub root {
    # Getter / setter for root node.
    my $self = shift;
    if (@_) {
       my $value = shift;
       if( defined $value && 
       ! $value->isa('Bio::Tree::NodeI') ) { 
       $self->warn("Trying to set the root node to $value which is not a Bio::Tree::NodeI");
       return $self->get_root_node;
       }
       $self->{'_rootnode'} = $value;
    }
    return $self->{'_rootnode'};
}
sub get_root_node { shift->root }
sub set_root_node { shift->root(@_) }

sub id{
   my ($self,$val) = @_;
   if( defined $val ) { 
       $self->{'_treeid'} = $val;
   }
   return $self->{'_treeid'};
}

sub score{
   my ($self,$val) = @_;
   if( defined $val ) { 
       $self->{'_score'} = $val;
   }
   return $self->{'_score'};
}

# safe tree clone that doesn't seg fault

=head2 clone()

 Title   : clone
 Alias   : _clone
 Usage   : $tree_copy = $tree->clone();
           $subtree_copy = $tree->clone($internal_node);
 Function: Safe tree clone that doesn't segfault
           (of Sendu)
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
        for my $node ($self->get_root_node->nodes) {
            $node->node_cleanup;
        }
    }
    $self->{'_rootnode'} = undef;
}

1;
