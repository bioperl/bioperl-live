# $Id$
#
# BioPerl module for Bio::Tree::Tree
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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

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


=head2 nodelete

 Title   : nodelete
 Usage   : $obj->nodelete($newval)
 Function: Get/Set Boolean whether or not to delete the underlying
           nodes when it goes out of scope.  By default this is false
           meaning trees are cleaned up.
 Returns : boolean
 Args    : on set, new boolean value


=cut

sub nodelete{
    my $self = shift;
    return $self->{'nodelete'} = shift if @_;
    return $self->{'nodelete'};
}

=head2 get_nodes

 Title   : get_nodes
 Usage   : my @nodes = $tree->get_nodes()
 Function: Return list of Tree::NodeI objects
 Returns : array of Tree::NodeI objects
 Args    : (named values) hash with one value 
           order => 'b|breadth' first order or 'd|depth' first order

=cut

sub get_nodes{
   my ($self, @args) = @_;
   
   my ($order, $sortby) = $self->_rearrange([qw(ORDER SORTBY)],@args);
   $order ||= 'depth';
   $sortby ||= 'none';
   my $node = $self->get_root_node || return;
   if ($order =~ m/^b|(breadth)$/oi) {
       my @children = ($node);
       for (@children) {
        push @children, $_->each_Descendent($sortby);
       }
       return @children;
   }

   if ($order =~ m/^d|(depth)$/oi) {
       # this is depth-first search I believe
       my @children = ($node,$node->get_all_Descendents($sortby));
       return @children;
   }
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
   return $self->{'_rootnode'};
}

=head2 set_root_node

 Title   : set_root_node
 Usage   : $tree->set_root_node($node)
 Function: Set the Root Node for the Tree
 Returns : Bio::Tree::NodeI
 Args    : Bio::Tree::NodeI

=cut

sub set_root_node{
   my $self = shift;
   if( @_ ) { 
       my $value = shift;
       if( defined $value && 
	   ! $value->isa('Bio::Tree::NodeI') ) { 
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
 Returns : integer
 Args    : none

=cut

sub total_branch_length {
   my ($self) = @_;
   my $sum = 0;
   if( defined $self->get_root_node ) {
       for ( $self->get_root_node->get_all_Descendents('none') ) {
	   $sum += $_->branch_length || 0;
       }
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

sub id{
   my ($self,$val) = @_;
   if( defined $val ) { 
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

sub score{
   my ($self,$val) = @_;
   if( defined $val ) { 
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


=cut


# -- private internal methods --

sub cleanup_tree {
    my $self = shift;
    unless( $self->nodelete ) {
        for my $node ($self->get_nodes(-order  => 'b', -sortby => 'none')) {
            #$node->ancestor(undef);
            #$node = undef;
            $node->node_cleanup;
            undef $node;
        }
    }
    $self->{'_rootnode'} = undef;
}

1;
