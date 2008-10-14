# $Id$
#
# BioPerl module for Bio::Tree::Statistics
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
previously where statistics from a Coalescent simulation.  Currently
it is empty because we have not added any Tree specific statistic
calculations to this module yet.  We welcome any contributions.


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


=cut


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
  Args       : Bio::Tree::TreeI object
               Bio::Tree::NodeI object within the tree, optional

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

The following methods produce desciptors of trait distribution among
leaf nodes within the trees. They require that a trait has to be set
for each leaf node. The tag methods of Bio::Tree::Node are used to
store them as key/value pairs. In this way, one tree can store more
than on trait.

Trees have method add_traits() to set trait values from a file.

=head2 ps

  Example    : ps($tree, $key, $node);
  Description: Calculates Parsimony Score (PS) from Fitch 1971
               parsimony algorithm for the subtree a defined
               by the (internal) node.
               Node defaults to the root.
  Returns    : integer, 1< PS < n, where n is number of branches
  Exceptions : leaf nodes have to have the trait defined
  Args       : Bio::Tree::TreeI object
               trait name string
               Bio::Tree::NodeI object within the tree, optional


Fitch, W.M., 1971. Toward deﬁning the course of evolution: minimal
change for a speciﬁc tree topology. Syst. Zool. 20, 406–416.

=cut

sub ps {
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
	$self->ps($tree, $key, $child);
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



=head2 ai

  Example    : ai($tree, $key, $node);
  Description: Calculates the Association Index (AI) of Whang et
               al. 2001 for the subtree defined by the (internal)
               node.  Node defaults to the root.
  Returns    : real
  Exceptions : leaf nodes have to have the trait defined
  Args       : Bio::Tree::TreeI object
               trait name string
               Bio::Tree::NodeI object within the tree, optional


Association index (AI) gives a more fine grained results than PS since
the result is a real number. ~0 <= AI.

Wang, T.H., Donaldson, Y.K., Brettle, R.P., Bell, J.E., Simmonds, P.,
  2001.  Identiﬁcation of shared populations of human immunodeﬁciency
  Virus Type 1 infecting microglia and tissue macrophages outside the
  central nervous system. J. Virol. 75 (23), 11686–11699.

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
  Args       : Bio::Tree::TreeI object
               trait name string
               Bio::Tree::NodeI object within the tree, optional


* Monophyletic Clade (MC) size statistics by Salemi at al 2005. It is
calculated for each trait value. 1<= MC <= nx, where nx is the
number of tips with value x:

   pick the internal node with maximim value for
      number of of tips with only trait x

MC was defined by Parker et al 2008.

Salemi, M., Lamers, S.L., Yu, S., de Oliveira, T., Fitch, W.M., McGrath, M.S.,
   2005. Phylodynamic analysis of Human Immunodeﬁciency Virus Type 1 in
   distinct brain compartments provides a model for the neuropathogenesis of
   AIDS. J. Virol. 79 (17), 11343–11352.

Parker, J., Rambaut A., Pybus O., 2008. Correlating viral phenotypes
   with phylogeny: Accounting for phylogenetic uncertainty Infection,
   Genetics and Evolution 8 (2008), 239–246.

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
