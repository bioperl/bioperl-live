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

Email jason@bioperl.org

=head1 CONTRIBUTORS

none so far

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
 Usage   : my $obj = new Bio::Tree::Statistics();
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


1;
