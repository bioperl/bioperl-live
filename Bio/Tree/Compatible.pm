#
# BioPerl module for Bio::Tree::Compatible
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gabriel Valiente <valiente@lsi.upc.edu>
#
# Copyright Gabriel Valiente
#
# You may distribute this module under the same terms as Perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Compatible - Testing compatibility of phylogenetic trees
with nested taxa.

=head1 SYNOPSIS

  use Bio::Tree::Compatible;
  use Bio::TreeIO;
  my $input = Bio::TreeIO->new('-format' => 'newick',
                               '-file'   => 'input.tre');
  my $t1 = $input->next_tree;
  my $t2 = $input->next_tree;

  my ($incompat, $ilabels, $inodes) = Bio::Tree::Compatible::is_compatible($t1,$t2);
  if ($incompat) {
    my %cluster1 = %{ Bio::Tree::Compatible::cluster_representation($t1) };
    my %cluster2 = %{ Bio::Tree::Compatible::cluster_representation($t2) };
    print "incompatible trees\n";
    if (scalar(@$ilabels)) {
      foreach my $label (@$ilabels) {
        my $node1 = $t1->find_node(-id => $label);
        my $node2 = $t2->find_node(-id => $label);
        my @c1 = sort @{ $cluster1{$node1} };
        my @c2 = sort @{ $cluster2{$node2} };
        print "label $label";
        print " cluster"; map { print " ",$_ } @c1;
        print " cluster"; map { print " ",$_ } @c2; print "\n";
      }
    }
    if (scalar(@$inodes)) {
      while (@$inodes) {
        my $node1 = shift @$inodes;
        my $node2 = shift @$inodes;
        my @c1 = sort @{ $cluster1{$node1} };
        my @c2 = sort @{ $cluster2{$node2} };
        print "cluster"; map { print " ",$_ } @c1;
        print " properly intersects cluster";
        map { print " ",$_ } @c2; print "\n";
      }
    }
  } else {
    print "compatible trees\n";
  }

=head1 DESCRIPTION

NB: This module has exclusively class methods that work on Bio::Tree::TreeI
objects. An instance of Bio::Tree::Compatible cannot itself represent a tree,
and so typically there is no need to create one.

Bio::Tree::Compatible is a Perl tool for testing compatibility of
phylogenetic trees with nested taxa represented as Bio::Tree::Tree
objects. It is based on a recent characterization of ancestral
compatibility of semi-labeled trees in terms of their cluster
representations.

A semi-labeled tree is a phylogenetic tree with some of its internal
nodes labeled, and it can represent a classification tree as well as a
phylogenetic tree with nested taxa, with labeled internal nodes
corresponding to taxa at a higher level of aggregation or nesting than
that of their descendents.

Two semi-labeled trees are compatible if their topological
restrictions to the common labels are such that for each node label,
the smallest clusters containing it in each of the trees coincide and,
furthermore, no cluster in one of the trees properly intersects a
cluster of the other tree.

Future extensions of Bio::Tree::Compatible include a
Bio::Tree::Supertree module for combining compatible phylogenetic
trees with nested taxa into a common supertree.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

=head1 SEE ALSO

=over

=item * Philip Daniel and Charles Semple. Supertree Algorithms for
Nested Taxa. In: Olaf R. P. Bininda-Emonds (ed.) Phylogenetic
Supertrees: Combining Information to Reveal the Tree of Life,
I<Computational Biology>, vol. 4, chap. 7, pp. 151-171. Kluwer (2004).

=item * Charles Semple, Philip Daniel, Wim Hordijk, Roderic
D. M. Page, and Mike Steel: Supertree Algorithms for Ancestral
Divergence Dates and Nested Taxa. Bioinformatics B<20>(15), 2355-2360
(2004).

=item * Merce Llabres, Jairo Rocha, Francesc Rossello, and Gabriel
Valiente: On the Ancestral Compatibility of Two Phylogenetic Trees
with Nested Taxa. J. Math. Biol. B<53>(3), 340-364 (2006).

=back

=head1 AUTHOR - Gabriel Valiente

Email valiente@lsi.upc.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::Tree::Compatible;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Set::Scalar;

use base qw(Bio::Root::Root);

=head2 postorder_traversal

 Title   : postorder_traversal
 Usage   : my @nodes = @{ $tree->postorder_traversal }
 Function: Return list of nodes in postorder
 Returns : reference to array of Bio::Tree::Node
 Args    : none

For example, the postorder traversal of the tree
C<(((A,B)C,D),(E,F,G));> is a reference to an array of nodes with
internal_id 0 through 9, because the Newick standard representation
for phylogenetic trees is based on a postorder traversal.

          +---A                    +---0
          |                        |
  +---+---C                +---4---2
  |   |   |                |   |   |
  |   |   +---B            |   |   +---1
  |   |                    |   |
  +   +-------D            9   +-------3
  |                        |
  |     +-----E            |     +-----5
  |     |                  |     |
  +-----+-----F            +-----8-----6
        |                        |
        +-----G                  +-----7

=cut

sub postorder_traversal {
  my($self) = @_;
  my @stack;
  my @queue;
  push @stack, $self->get_root_node;
  while (@stack) {
    my $node = pop @stack;
    push @queue, $node;
    foreach my $child ($node->each_Descendent(-sortby => 'internal_id')) {
      push @stack, $child;
    }
  }
  my @postorder = reverse @queue;
  return \@postorder;
}

=head2 cluster_representation

 Title   : cluster_representation
 Usage   : my %cluster = %{ $tree->cluster_representation }
 Function: Compute the cluster representation of a tree
 Returns : reference to hash of array of string indexed by
           Bio::Tree::Node
 Args    : none

For example, the cluster representation of the tree
C<(((A,B)C,D),(E,F,G));> is a reference to a hash associating an array
of string (descendent labels) to each node, as follows:

  0 --> [A]
  1 --> [B]
  2 --> [A,B,C]
  3 --> [D]
  4 --> [A,B,C,D]
  5 --> [E]
  6 --> [F]
  7 --> [G]
  8 --> [E,F,G]
  9 --> [A,B,C,D,E,F,G]

=cut

sub cluster_representation {
  my ($tree) = @_;
  my %cluster;
  my @postorder = @{ postorder_traversal($tree) };
  foreach my $node ( @postorder ) {
    my @labeled = map { $_->id } grep { $_->id } $node->get_Descendents;
    push @labeled, $node->id if $node->id;
    $cluster{$node} = \@labeled;
  }
  return \%cluster;
}

=head2 common_labels

 Title   : common_labels
 Usage   : my $labels = $tree1->common_labels($tree2);
 Function: Return set of common node labels
 Returns : Set::Scalar
 Args    : Bio::Tree::Tree

For example, the common labels of the tree C<(((A,B)C,D),(E,F,G));>
and the tree C<((A,B)H,E,(J,(K)G)I);> are: C<[A,B,E,G]>.

          +---A                 +---A
          |                     |
  +---+---C             +-------H
  |   |   |             |       |
  |   |   +---B         |       +---B
  |   |                 |
  +   +-------D         +-----------E
  |                     |
  |     +-----E         |   +-------J
  |     |               |   |
  +-----+-----F         +---I
        |                   |
        +-----G             +---G---K

=cut

sub common_labels {
  my($self,$arg) = @_;
  my @labels1 = map { $_->id } grep { $_->id } $self->get_nodes;
  my $common = Set::Scalar->new( @labels1 );
  my @labels2 = map { $_->id } grep { $_->id } $arg->get_nodes;
  my $temp = Set::Scalar->new( @labels2 );
  return $common->intersection($temp);
}

=head2 topological_restriction

 Title   : topological_restriction
 Usage   : $tree->topological_restriction($labels)
 Function: Compute the topological restriction of a tree to a subset
           of node labels
 Returns : Bio::Tree::Tree
 Args    : Set::Scalar

For example, the topological restrictions of each of the trees
C<(((A,B)C,D),(E,F,G));> and C<((A,B)H,E,(J,(K)G)I);> to the labels
C<[A,B,E,G]> are as follows:

          +---A             +---A
          |                 |
  +---+---+             +---+
  |       |             |   |
  |       +---B         |   +---B
  +                     |
  |       +---E         +-------E
  |       |             |
  +-------+             +---+---G
          |
          +---G

=cut

sub topological_restriction {
  my ($tree, $labels) = @_;
  for my $node ( @{ postorder_traversal($tree) } ) {
    unless (ref($node)) { # skip $node if already removed
      my @cluster = map { $_->id } grep { $_->id } $node->get_Descendents;
      push @cluster, $node->id if $node->id;
      my $cluster = Set::Scalar->new(@cluster);
      if ($cluster->is_disjoint($labels)) {
        $tree->remove_Node($node);
      } else {
        if ($node->id and not $labels->has($node->id)) {
          $node->{'_id'} = undef;
        }
      }
    }
  }
}

=head2 is_compatible

 Title   : is_compatible
 Usage   : $tree1->is_compatible($tree2)
 Function: Test compatibility of two trees
 Returns : boolean
 Args    : Bio::Tree::Tree

For example, the topological restrictions of the trees
C<(((A,B)C,D),(E,F,G));> and C<((A,B)H,E,(J,(K)G)I);> to their common
labels, C<[A,B,E,G]>, are compatible. The respective cluster
representations are as follows:

  [A]                  [A]
  [B]                  [B]
  [E]                  [E]
  [G]                  [G]
  [A,B]                [A,B]
  [E,G]                [A,B,E,G]
  [A,B,E,G]

As a second example, the trees C<(A,B);> and C<((B)A);> are
incompatible. Their respective cluster representations are as follows:

  [A]                  [B]
  [B]                  [A,B]
  [A,B]

The reason is, the smallest cluster containing label C<A> is C<[A]> in
the first tree but C<[A,B]> in the second tree.

 +---A         A---B
 |
 +
 |
 +---B

As a second example, the trees C<(((B,A),C),D);> and C<((A,(D,B)),C);>
are also incompatible. Their respective cluster representations are as
follows:

  [A]                  [A]
  [B]                  [B]
  [C]                  [C]
  [D]                  [D]
  [A,B]                [B,D]
  [A,B,C]              [A,B,D]
  [A,B,C,D]            [A,B,C,D]

The reason is, cluster C<[A,B]> properly intersects cluster
C<[B,D]>. There are further incompatibilities between these trees:
C<[A,B,C]> properly intersects both C<[B,D]> and C<[A,B,D]>.

          +---B             +-------A
          |                 |
      +---+             +---+   +---D
      |   |             |   |   |
  +---+   +---A         |   +---+
  |   |                 +       |
  +   +-------C         |       +---B
  |                     |
  +-----------D         +-----------C

=cut

sub is_compatible {
  my ($tree1, $tree2) = @_;
  my $common = $tree1->Bio::Tree::Compatible::common_labels($tree2);
  $tree1->Bio::Tree::Compatible::topological_restriction($common);
  $tree2->Bio::Tree::Compatible::topological_restriction($common);
  my @postorder1 = @{ postorder_traversal($tree1) };
  my @postorder2 = @{ postorder_traversal($tree2) };
  my %cluster1 = %{ cluster_representation($tree1) };
  my %cluster2 = %{ cluster_representation($tree2) };
  my $incompat = 0; # false
  my @labels;
  foreach my $label ( $common->elements ) {
    my $node1 = $tree1->find_node(-id => $label);
    my @labels1 = @{ $cluster1{$node1} };
    my $cluster1 = Set::Scalar->new(@labels1);
    my $node2 = $tree2->find_node(-id => $label);
    my @labels2 = @{ $cluster2{$node2} };
    my $cluster2 = Set::Scalar->new(@labels2);
    unless ( $cluster1->is_equal($cluster2) ) {
      $incompat = 1; # true
      push @labels, $label;
    }
  }
  my @nodes;
  foreach my $node1 ( @postorder1 ) {
    my @labels1 = @{ $cluster1{$node1} };
    my $cluster1 = Set::Scalar->new(@labels1);
    foreach my $node2 ( @postorder2 ) {
      my @labels2 = @{$cluster2{$node2} };
      my $cluster2 = Set::Scalar->new(@labels2);
      if ($cluster1->is_properly_intersecting($cluster2)) {
	$incompat = 1; # true
	push @nodes, $node1, $node2;
      }
    }
  }
  return ($incompat, \@labels, \@nodes);
}

1;
