#
# BioPerl module for Bio::ClusterI
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

# see the implementations of this interface for details but
# basically

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org


=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::ClusterI;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

=head1 Implementation Specific Functions

These functions are the ones that a specific implementation must
define.

=head2 new

 Title   : new
 Usage   : Bio::ClusterI->new(@args)
 Function: Create new Bio::ClusterI object
 Returns : a L<Bio::ClusterI> object
 Args    : see below

  Argument        Description
  --------        -----------
  -description    the consensus description of the family
  -members        the array of objects belonging to the family

=cut



=head2 description

 Title   : description
 Usage   : Bio::ClusterI->description("POLYUBIQUITIN")
 Function: get/set for the consensus description of the cluster
 Returns : the description string 
 Args    : Optional the description string 

=cut

sub description{
  my ($self,$desc) = @_;
  if($desc){
    $self->{'_description'} = $desc;
  }
  return $self->{'_description'};
}

=head2 members

 Title   : members
 Usage   : Bio::ClusterI->members(($seq1, $seq2));
 Function: get/set for the members of the family 
 Returns : the array of members 
 Args    : the array of members 

=cut

sub members {
  my ($self,$mem) = @_;
  if($mem){
      $self->{'_members'} = ();
      if(ref($mem) eq "ARRAY"){
          push @{$self->{'_members'}}, @{$mem};
      }
      else {
        push @{$self->{'_members'}}, $mem;
      }
  }
  return @{$self->{'_members'}};

}

=head2 add_members

 Title   : add_members
 Usage   : Bio::ClusterI->add_member([$seq1,$seq1]);
 Function: add members to a family
 Returns : 
 Args    : the arrayref of members to add or a single member

=cut

sub add_members{
  my ($self,$mem) = @_;
  if($mem){
      if(ref($mem) eq "ARRAY"){
        push @{$self->{'_members'}},@{$mem};
      }
      else {
        push @{$self->{'_members'}},$mem;
      }
  }

  return @{$self->{'_members'}};
}

*add_member = \&add_members;




=head2 flush_members

 Title   : flush_members
 Usage   : Bio::ClusterI->flush_members();
 Function: remove all members from a family 
 Returns :
 Args    : 

=cut

sub flush_members{
    my ($self) =  @_;
    $self->{'_members'} = [];
    
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
  my ($self) = @_;

  return scalar(@{$self->{'_members'}});

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
  my ($self,$score) = @_;
  if($score){
    $self->{'_cluster_score'} = $score;
  }
  return $self->{'_cluster_score'};
}

############################
#Abstract Methods
############################

=head2 get_members

 Title   : get_members
 Usage   : Bio::ClusterI->get_members(($seq1, $seq2));
 Function: retrieve the members of the family by some criteria, for
           example :
           $cluster->get_members(-species => 'homo sapiens'); 
 Returns : the array of members
 Args    : 

=cut

sub get_members {
  shift->throw_not_implemented();
}

=head2 remove_members

 Title   : remove_members
 Usage   : Bio::ClusterI->remove_members(("ENSP000000001","ENSP000000002"));
 Function: remove members from a family
 Returns :
 Args    : the display_ids of the members 

=cut

sub remove_members{
  shift->throw_not_implemented();
}
1;
