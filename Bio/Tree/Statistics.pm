#
# BioPerl module for Bio::Tree::Statistics
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

Bio::Tree::Statistics - Calculate certain statistics for a Tree

=head1 SYNOPSIS

  use Bio::Tree::Statistics;

=head1 DESCRIPTION

This should be where Tree statistics are calculated.  It was
previously where statistics from a Coalescent simulation.

It now contains several methods for calculating L<Tree-Trait
statistics>. 

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

Email jason AT bioperl.org

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::Statistics;
use strict;


use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tree::Statistics->new();
 Function: Builds a new Bio::Tree::Statistics object 
 Returns : Bio::Tree::Statistics
 Args    :

=head2 assess_bootstrap

 Title   : assess_bootstrap
 Usage   : my $tree_with_bs = $stats->assess_bootstrap(\@bs_trees);
 Function: Calculates the bootstrap for internal nodes based on
 Returns : L<Bio::Tree::TreeI>
 Args    : Arrayref of L<Bio::Tree::TreeI>s

=cut

sub assess_bootstrap{
   my ($self,$bs_trees,$guide_tree) = @_;
   my @consensus;

   # internal nodes are defined by their children

   my (%lookup,%internal);
   my $i = 0;
   for my $tree ( $guide_tree, @$bs_trees ) {
       # Do this as a top down approach, can probably be
       # improved by caching internal node states, but not going
       # to worry about it right now.

       my @allnodes = $tree->get_nodes;
       my @internalnodes = grep { ! $_->is_Leaf } @allnodes;
       for my $node ( @internalnodes ) {
           my @tips = sort map { $_->id } 
                      grep { $_->is_Leaf() } $node->get_all_Descendents;
           my $id = "(".join(",", @tips).")";
           if( $i == 0 ) {
               $internal{$id} = $node->internal_id;
           } else { 
               $lookup{$id}++;
           }
       }
       $i++;
   }
   my @save;
   for my $l ( keys %lookup ) {
       if( defined $internal{$l} ) {#&& $lookup{$l} > $min_seen ) {
           my $intnode = $guide_tree->find_node(-internal_id => $internal{$l});
           $intnode->bootstrap(sprintf("%d",100 * $lookup{$l} / $i));
       }
   }
   return $guide_tree;
}


=head2 cherries

  Example    : cherries($tree, $node);
  Description: Count number of paired leaf nodes
               in a binary tree
  Returns    : integer
  Exceptions : 
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Commonly used statistics assume a binary tree, but this methods
returns a value even for trees with polytomies.

=cut

sub cherries ($;$) {
    my $self = shift;
    my $tree = shift;
    my $node = shift || $tree->get_root_node;

    my $cherries = 0;
    my @descs = $node->each_Descendent;

    if ($descs[0]->is_Leaf and $descs[1]->is_Leaf) {
        if ($descs[3]) { #polytomy at leaf level
            $cherries = 0;
        } else {
            $cherries = 1;
        }
    } else {
        # recurse
        foreach my $desc (@descs) {
            $cherries += $self->cherries($tree, $desc);
        }
    }
    return $cherries;
}


=head2 Tree-Trait statistics

The following methods produce descriptors of trait distribution among
leaf nodes within the trees. They require that a trait has been set
for each leaf node. The tag methods of Bio::Tree::Node are used to
store them as key/value pairs. In this way, one tree can store more
than one trait.

Trees have method add_traits() to set trait values from a file. See the
add_trait() method in L<Bio::Tree::TreeFunctionsI>.

=head2 fitch

  Example    : fitch($tree, $key, $node);
  Description: Calculates Parsimony Score (PS) and internal trait
               values using the Fitch 1971 parsimony algorithm for
               the subtree a defined by the (internal) node.
               Node defaults to the root.
  Returns    : true on success
  Exceptions : leaf nodes have to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. trait name string
               3. Bio::Tree::NodeI object within the tree, optional

Runs first L<fitch_up> that calculates parsimony scores and then
L<fitch_down> that should resolve most of the trait/character state
ambiguities.

Fitch, W.M., 1971. Toward defining the course of evolution: minimal
change for a specific tree topology. Syst. Zool. 20, 406-416.

You can access calculated parsimony values using:

  $score = $node->->get_tag_values('ps_score');

and the trait value with:

  $traitvalue = $node->->get_tag_values('ps_trait'); # only the first
  @traitvalues = $node->->get_tag_values('ps_trait');

Note that there can be more that one trait value, especially for the
root node.

=cut

sub fitch {
    my $self = shift;
    my $tree = shift;
    my $key = shift || $self->throw("Trait name is needed");
    my $node = shift || $tree->get_root_node;

    $self->fitch_up($tree, $key, $node);
    $self->fitch_down($tree, $node);
}


=head2 ps

  Example    : ps($tree, $key, $node);
  Description: Calculates Parsimony Score (PS) from Fitch 1971
               parsimony algorithm for the subtree as defined
               by the (internal) node.
               Node defaults to the root.
  Returns    : integer, 1 < PS < n, where n is number of branches
  Exceptions : leaf nodes have to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. trait name string
               3. Bio::Tree::NodeI object within the tree, optional

This is the first half of the Fitch algorithm that is enough for
calculating the resolved parsimony values. The trait/chararacter
states are commonly left in ambiguous state. To resolve them, run
L<fitch_down>.

=cut

sub ps { shift->fitch_up(@_) }


=head2 fitch_up

  Example    : fitch_up($tree, $key, $node);
  Description: Calculates Parsimony Score (PS) from the Fitch 1971
               parsimony algorithm for the subtree as defined
               by the (internal) node.
               Node defaults to the root.
  Returns    : integer, 1< PS < n, where n is number of branches
  Exceptions : leaf nodes have to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. trait name string
               3. Bio::Tree::NodeI object within the tree, optional

This is a more generic name for L<ps> and indicates that it performs
the first bottom-up tree traversal that calculates the parsimony score
but usually leaves trait/character states ambiguous. If you are
interested in internal trait states, running L<fitch_down> should
resolve most of the ambiguities.

=cut

sub fitch_up {
    my $self = shift;
    my $tree = shift;
    my $key = shift || $self->throw("Trait name is needed");
    my $node = shift || $tree->get_root_node;

    if ($node->is_Leaf) {
        $self->throw ("ERROR: ". $node->internal_id. " needs a value for trait $key")
            unless $node->has_tag($key);
        $node->set_tag_value('ps_trait', $node->get_tag_values($key) );
        $node->set_tag_value('ps_score', 0 );
        return; # end of recursion
    }

    foreach my $child ($node->each_Descendent) {
        $self->fitch_up($tree, $key, $child);
    }

    my %intersection;
    my %union;
    my $score;

    foreach my $child ($node->each_Descendent) {
        foreach my $trait ($child->get_tag_values('ps_trait') ) {
            $intersection{$trait}++ if $union{$trait};
            $union{$trait}++;
        }
        $score += $child->get_tag_values('ps_score');
    }

    if (keys %intersection) {
        $node->set_tag_value('ps_trait', keys %intersection);
        $node->set_tag_value('ps_score', $score);
    } else {
        $node->set_tag_value('ps_trait', keys %union);
        $node->set_tag_value('ps_score', $score+1);
    }

    if ($self->verbose) {
        print "-- node --------------------------\n";
        print "iID: ", $node->internal_id, " (", $node->id, ")\n";
        print "Trait: ", join (', ', $node->get_tag_values('ps_trait') ), "\n";
        print "length :", scalar($node->get_tag_values('ps_score')) , "\n";
    }
    return scalar $node->get_tag_values('ps_score');
}


=head2 fitch_down

  Example    : fitch_down($tree, $node);
  Description: Runs the second pass from Fitch 1971
               parsimony algorithm to resolve ambiguous
               trait states left by first pass.
               by the (internal) node.
               Node defaults to the root.
  Returns    : true
  Exceptions : dies unless the trait is defined in all nodes
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Before running this method you should have ran L<fitch_up> (alias to
L<ps> ). Note that it is not guaranteed that all states are completely
resolved.

=cut

sub fitch_down {
    my $self = shift;
    my $tree = shift;
    my $node = shift || $tree->get_root_node;

    my $key = 'ps_trait';
    $self->throw ("ERROR: ". $node->internal_id. " needs a value for $key")
        unless $node->has_tag($key);

    my $nodev;
    foreach my $trait ($node->get_tag_values($key) ) {
        $nodev->{$trait}++;
    }

    foreach my $child ($node->each_Descendent) {
        next if $child->is_Leaf;  # end of recursion

        my $intersection;
        foreach my $trait ($child->get_tag_values($key) ) {
            $intersection->{$trait}++ if $nodev->{$trait};
        }

        $self->fitch_down($tree, $child);
        $child->set_tag_value($key, keys %$intersection);
    }
    return 1;  # success
}


=head2 persistence

  Example    : persistence($tree, $node);
  Description: Calculates the persistence
               for node in the subtree defined by the (internal)
               node.  Node defaults to the root.
  Returns    : int, number of generations trait value has to remain same
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Persistence measures the stability that the trait value has in a
tree. It expresses the number of generations the trait value remains
the same. All the decendants of the root in the same generation have
to share the same value.

Depends on Fitch's parsimony score (PS).

=cut

sub _persistence {
    my $self = shift;
    my $tree = shift;
    my $node = shift;
    my $value = shift || $self->throw("Value is needed");


    my $key  = 'ps_trait';

    $self->throw("Node is needed") unless $node->isa('Bio::Tree::NodeI');

    return 0 unless $node->get_tag_values($key) eq $value; # wrong value
    return 1 if $node->is_Leaf; # end of recursion

    my $persistence = 10000000; # an arbitrarily large number
    foreach my $child ($node->each_Descendent) {
        my $pers = $self->_persistence($tree, $child, $value);
        $persistence = $pers if $pers < $persistence;
    }
    return $persistence + 1;
}

sub persistence {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;
    $self->throw("Node is needed") unless $node->isa('Bio::Tree::NodeI');

    my $key  = 'ps_trait';
    my $value = $node->get_tag_values($key);

    #calculate
    my $persistence =  $self->_persistence($tree, $node, $value);
    $node->set_tag_value('persistance', $persistence);
    return $persistence;
}


=head2 count_subclusters

  Example    : count_clusters($tree, $node);
  Description: Calculates the number of sub-clusters
               in the subtree defined by the (internal)
               node.  Node defaults to the root.
  Returns    : int, count
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Depends on Fitch's parsimony score (PS).

=cut

sub _count_subclusters {
    my $self = shift;
    my $tree = shift;
    my $node = shift;
    my $value = shift || $self->throw("Value is needed");

    my $key  = 'ps_trait';

    $self->throw ("ERROR: ". $node->internal_id. " needs a value for trait $key")
        unless $node->has_tag($key);

    if ($node->get_tag_values($key) eq $value) {
        if ($node->get_tag_values('ps_score') == 0) {
            return 0;
        } else {
            my $count = 0;
            foreach my $child ($node->each_Descendent) {
                $count += $self->_count_subclusters($tree, $child, $value);
            }
            return $count;
        }
    }
    return 1;
}

sub count_subclusters {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;
    $self->throw("Node is needed") unless $node->isa('Bio::Tree::NodeI');

    my $key  = 'ps_trait';
    my $value = $node->get_tag_values($key);

    return $self->_count_subclusters($tree, $node, $value);
}


=head2 count_leaves

  Example    : count_leaves($tree, $node);
  Description: Calculates the number of leaves with same trait
               value as root in the subtree defined by the (internal)
               node.  Requires an unbroken line of identical trait values.
               Node defaults to the root.
  Returns    : int, number of leaves with this trait value
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Depends on Fitch's parsimony score (PS).

=cut

sub _count_leaves {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;
    my $value = shift;

    my $key  = 'ps_trait';

    $self->throw ("ERROR: ". $node->internal_id. " needs a value for trait $key")
        unless $node->has_tag($key);

    if ($node->get_tag_values($key) eq $value) {
        #print $node->id, ": ", $node->get_tag_values($key), "\n";
        return 1 if $node->is_Leaf; # end of recursion

            my $count = 0;
            foreach my $child ($node->each_Descendent) {
                $count += $self->_count_leaves($tree, $child, $value);
            }
            return $count;
    }
    return 0;
}

sub count_leaves {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;
    $self->throw("Node is needed") unless $node->isa('Bio::Tree::NodeI');

    my $key  = 'ps_trait';
    my $value = $node->get_tag_values($key);

    return $self->_count_leaves($tree, $node, $value);
}


=head2 phylotype_length

  Example    : phylotype_length($tree, $node);
  Description: Sums up the branch lengths within phylotype
               exluding the subclusters where the trait values
               are different
  Returns    : float, length
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Depends on Fitch's parsimony score (PS).

=cut

sub _phylotype_length {
    my $self = shift;
    my $tree = shift;
    my $node = shift;
    my $value = shift;

    my $key  = 'ps_trait';

    $self->throw ("ERROR: ". $node->internal_id. " needs a value for trait $key")
        unless $node->has_tag($key);

    return 0 if $node->get_tag_values($key) ne $value;
    return $node->branch_length if $node->is_Leaf; # end of recursion

    my $length = 0;
    foreach my $child ($node->each_Descendent) {
        my $sub_len = $self->_phylotype_length($tree, $child, $value);
        $length += $sub_len;
        $length += $child->branch_length if not $child->is_Leaf and $sub_len;
    }
    return $length;
}

sub phylotype_length {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;

    my $key  = 'ps_trait';
    my $value = $node->get_tag_values($key);

    return $self->_phylotype_length($tree, $node, $value);
}


=head2 sum_of_leaf_distances

  Example    : sum_of_leaf_distances($tree, $node);
  Description: Sums up the branch lengths from root to leaf
               exluding the subclusters where the trait values
               are different
  Returns    : float, length
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Depends on Fitch's parsimony score (PS).

=cut

sub _sum_of_leaf_distances {
    my $self = shift;
    my $tree = shift;
    my $node = shift;
    my $value = shift;

    my $key  = 'ps_trait';

    $self->throw ("ERROR: ". $node->internal_id. " needs a value for trait $key")
        unless $node->has_tag($key);
    return 0 if $node->get_tag_values($key) ne $value;
    #return $node->branch_length if $node->is_Leaf; # end of recursion
    return 0 if $node->is_Leaf; # end of recursion

    my $length = 0;
    foreach my $child ($node->each_Descendent) {
        $length += $self->_count_leaves($tree, $child, $value) * $child->branch_length +
        $self->_sum_of_leaf_distances($tree, $child, $value);
    }
    return $length;
}

sub sum_of_leaf_distances {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;

    my $key  = 'ps_trait';
    my $value = $node->get_tag_values($key);

    return $self->_sum_of_leaf_distances($tree, $node, $value);
}


=head2 genetic_diversity

  Example    : genetic_diversity($tree, $node);
  Description: Diversity is the sum of root to leaf distances
               within the phylotype normalised by number of leaf
               nodes
  Returns    : float, value of genetic diversity
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Depends on Fitch's parsimony score (PS).

=cut

sub genetic_diversity {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;

    return $self->sum_of_leaf_distances($tree, $node) /
        $self->count_leaves($tree, $node);
}


=head2 statratio

  Example    : statratio($tree, $node);
  Description: Ratio of the stem length and the genetic diversity of the
               phylotype L<genetic_diversity>
  Returns    : float, separation score
  Exceptions : all the  nodes need to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. Bio::Tree::NodeI object within the tree, optional

Statratio gives a measure of separation and variability within the phylotype.
Larger values identify more rapidly evolving and recent phylotypes.

Depends on Fitch's parsimony score (PS).

=cut

sub statratio {
    my $self = shift;
    my $tree = shift;
    my $node = shift  || $tree->get_root_node;

    my $div = $self->genetic_diversity($tree, $node);
    return 0 if $div == 0;
    return $node->branch_length / $div;

}


=head2 ai

  Example    : ai($tree, $key, $node);
  Description: Calculates the Association Index (AI) of Whang et
               al. 2001 for the subtree defined by the (internal)
               node.  Node defaults to the root.
  Returns    : real
  Exceptions : leaf nodes have to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. trait name string
               3. Bio::Tree::NodeI object within the tree, optional

  Association index (AI) gives a more fine grained results than PS since
  the result is a real number. ~0 E<lt>= AI.

  Wang, T.H., Donaldson, Y.K., Brettle, R.P., Bell, J.E., Simmonds, P.,
  2001.  Identification of shared populations of human immunodeficiency
  Virus Type 1 infecting microglia and tissue macrophages outside the
  central nervous system. J. Virol. 75 (23), 11686-11699.

=cut

sub _node_ai {
    my $self = shift;
    my $node = shift;
    my $key = shift;

    my $traits;
    my $leaf_count = 0;
    for my $desc ( $node->get_all_Descendents ) {
        next unless $desc->is_Leaf;
        $leaf_count++;
        $self->throw ("Node ". $desc->id. " needs a value for trait [$key]")
            unless $desc->has_tag($key);
        my $trait = $desc->get_tag_values($key);
        $traits->{$trait}++;
    }
    my $most_common = 0;
    foreach ( keys %$traits) {
        $most_common = $traits->{$_} if $traits->{$_} > $most_common;
    }
    return sprintf "%1.6f", (1 - ($most_common/$leaf_count) ) / (2**($leaf_count-1) );
}

sub ai {
    my $self = shift;
    my $tree = shift;
    my $key = shift || $self->throw("Trait name is needed");
    my $start_node = shift || $tree->get_root_node;
    return unless $start_node;

    my $sum = 0;
    for my $node ( $start_node->get_all_Descendents ) {
        next if $node->is_Leaf;
        $sum += $self->_node_ai($node, $key);
    }
    return $sum;
}


=head2 mc

  Example    : mc($tree, $key, $node);
  Description: Calculates the Monophyletic Clade (MC) size statistics
               for the subtree a defined by the (internal) node.
               Node defaults to the root;
  Returns    : hashref with trait values as keys
  Exceptions : leaf nodes have to have the trait defined
  Args       : 1. Bio::Tree::TreeI object
               2. trait name string
               3. Bio::Tree::NodeI object within the tree, optional

  Monophyletic Clade (MC) size statistics by Salemi at al 2005. It is
  calculated for each trait value. 1 E<lt>= MC E<lt>= nx, where nx is the
  number of tips with value x:

   pick the internal node with maximim value for
      number of of tips with only trait x

  MC was defined by Parker et al 2008.

  Salemi, M., Lamers, S.L., Yu, S., de Oliveira, T., Fitch, W.M., McGrath, M.S.,
   2005. Phylodynamic analysis of Human Immunodeficiency Virus Type 1 in
   distinct brain compartments provides a model for the neuropathogenesis of
   AIDS. J. Virol. 79 (17), 11343-11352.

  Parker, J., Rambaut A., Pybus O., 2008. Correlating viral phenotypes
   with phylogeny: Accounting for phylogenetic uncertainty Infection,
   Genetics and Evolution 8 (2008), 239-246.

=cut

sub _node_mc  {
    my $self = shift;
    my $node = shift;
    my $key = shift;

    my $traits;
    my $leaf_count = 0;
    for my $node2 ( $node->get_all_Descendents ) {
        next unless $node2->is_Leaf;
        $leaf_count++;
        my $trait = $node2->get_tag_values($key);
        $traits->{$trait}++;
    }
    return $traits;
}

sub mc {
    my $self = shift;
    my $tree = shift;
    my $key = shift || die "Trait name is needed";
    my $start_node = shift || $tree->get_root_node;
    return unless $start_node;

    my $sum; # hashref, keys are trait values
    my $keys; # hashref, keys are trait values
    foreach my $node ( $start_node->get_all_Descendents ) {
        next if $node->is_Leaf;
        my $traits = $self->_node_mc($node, $key);
        if (scalar keys %$traits == 1) {
            my ($value) = keys %$traits;
            no warnings;
            $sum->{$value} = $traits->{$value}
                if $sum->{$value} < $traits->{$value};
        } else {
            map { $keys->{$_} = 1 } keys %$traits;
        }
    }
    # check for cases where there are no clusters
    foreach my $value (keys %$keys) {
        $sum->{$value} = 1 unless defined $sum->{$value};
    }
    return $sum;
}


1;
