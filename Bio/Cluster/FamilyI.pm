#
# BioPerl module for Bio::Cluster::FamilyI
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Cluster::FamilyI - Family Interface

=head1 SYNOPSIS

# see the implementations of this interface for details but
# basically

    my $cluster= $cluster->new(-description=>"POLYUBIQUITIN",
                               -members    =>[$seq1,$seq2]);
    my @members = $cluster->get_members();
    my @sub_members = $cluster->get_members(-species=>"homo sapiens");



=head1 DESCRIPTION

This interface if for a Family object representing a family of 
biological objects. A generic implementation for this may be
found a Bio::Cluster::Family.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org


=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::Cluster::FamilyI;
use vars qw(@ISA);
use strict;

use Bio::ClusterI;

@ISA = qw(Bio::ClusterI);

=head2 new

 Title   : new
 Usage   : Bio::Cluster::FamilyI->new(@args)
 Function: Create new Bio::Cluster::FamilyI object
 Returns : a L<Bio::Cluster::FamilyI> object
 Args    : See below

 Arguments          Description
 ---------          -----------
 -family_id         the name of the family
 -description       the consensus description of the family
 -annotation_score  the confidence by which the consensus description is 
                    representative of the family
 -members           the members belonging to the family
 -alignment         the multiple alignment of the members

=cut



=head2 family_id

 Title   : family_id
 Usage   : Bio::Cluster::FamilyI->family_id("znfp");
 Function: get/set for the family id 
 Returns : the family id 
 Args    : the family id

=cut

sub family_id{
  my ($self,$id) = @_;
  if($id){
    $self->{'_family_id'} = $id;
  }
  return $self->{'_family_id'};
}

=head2 family_score

 Title   : family_score
 Usage   : Bio::Cluster::FamilyI->family_score(95);
 Function: get/set for the score of algorithm used to generate
           the family if present
 Returns : the score
 Args    : the score

=cut

sub family_score {
  my ($self,$score) = @_;
  if ($score){
    $self->cluster_score($score);
  }
  return $self->cluster_score;
}


=head1 Methods inherited from ClusterI

=head2 members

 Title   : members
 Usage   : Bio::Cluster::FamilyI->members(($seq1, $seq2));
 Function: get/set for the members of the family
 Returns : the array of members
 Args    : the array of members

=cut

=head2 description

 Title   : description
 Usage   : Bio::Cluster::FamilyI->description("Zinc Finger Protein");
 Function: get/set for the description of the family
 Returns : the description 
 Args    : the description

=cut


=head2 size

 Title   : size
 Usage   : Bio::Cluster::FamilyI->size();
 Function: get/set for the description of the family
 Returns : size 
 Args    : 

=cut


=head2 flush_members

 Title   : flush_members
 Usage   : Bio::Cluster::FamilyI->flush_members();
 Function: remove all members from a family
 Returns :
 Args    :

=cut

=head2 add_members

 Title   : add_members
 Usage   : Bio::Cluster::FamilyI->add_member(($seq1,$seq1));
 Function: add members to a family
 Returns :
 Args    : the array of members to add

=cut

=head2 cluster_score

 Title   : cluster_score
 Usage   : $cluster ->cluster_score(100);
 Function: get/set for cluster_score which
           represent the score in which the clustering
           algorithm assigns to this cluster.
 Returns : a number

=cut

1;
