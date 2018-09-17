#
# BioPerl module for Bio::ClusterI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::ClusterI - Cluster Interface 

=head1 SYNOPSIS

    # see the implementations of this interface for details

    my $cluster= $cluster->new(-description=>"POLYUBIQUITIN",
                               -members    =>[$seq1,$seq2]);
    my @members = $cluster->get_members();
    my @sub_members = $cluster->get_members(-species=>"homo sapiens");


=head1 DESCRIPTION

This interface is the basic structure for a cluster of bioperl objects.
In this case it is up to the implementer to check arguments
and initialize whatever new object the implementing class is designed for.

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

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::ClusterI;
use strict;


use base qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 new

  We don't mandate but encourage implementors to support at least the
  following named parameters upon object initialization.

  Argument        Description
  --------        -----------
  -display_id     the display ID or name for the cluster
  -description    the consensus description or name of the cluster
  -members        the array of objects belonging to the family

=cut

=head2 display_id

 Title   : display_id
 Usage   : 
 Function: Get the display name or identifier for the cluster
 Returns : a string
 Args    : 

=cut

sub display_id{
    shift->throw_not_implemented();
}


=head2 description

 Title   : description
 Usage   : Bio::ClusterI->description("POLYUBIQUITIN")
 Function: get/set for the consensus description of the cluster
 Returns : the description string 
 Args    : Optional the description string 

=cut

sub description{
    shift->throw_not_implemented();
}

=head2 size

 Title   : size
 Usage   : Bio::ClusterI->size();
 Function: get/set for the size of the family, 
           calculated from the number of members
 Returns : the size of the family 
 Args    : 

=cut

sub size {
    shift->throw_not_implemented();
}

=head2 cluster_score

 Title   : cluster_score
 Usage   : $cluster ->cluster_score(100);
 Function: get/set for cluster_score which
           represent the score in which the clustering
           algorithm assigns to this cluster.
 Returns : a number

=cut

sub cluster_score{
    shift->throw_not_implemented();
}

=head2 get_members

 Title   : get_members
 Usage   : Bio::ClusterI->get_members(($seq1, $seq2));
 Function: retrieve the members of the family by some criteria, for
           example :
           $cluster->get_members(-species => 'homo sapiens'); 

           Will return all members if no criteria are provided.

 Returns : the array of members
 Args    : 

=cut

sub get_members {
    shift->throw_not_implemented();
}

1;
