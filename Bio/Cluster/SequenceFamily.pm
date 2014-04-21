#
# BioPerl module for Bio::Cluster::SequenceFamily
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

Bio::Cluster::SequenceFamily - Sequence Family object

=head1 SYNOPSIS

  use Bio::SeqIO;
  use Bio::Cluster::SequenceFamily;
  use File::Spec;

  my $file =  File::Spec->catfile('t','data','swiss.dat');
  my $seqio= Bio::SeqIO->new(-format => 'swiss',
                            -file => $file);
  my @mem;
  while(my $seq = $seqio->next_seq){
    push @mem, $seq;
  }

  #create the family
  my $family = Bio::Cluster::SequenceFamily->new(
          -family_id=>"Family_1",
          -description=>"Family Description Here",
          -annotation_score=>"100",
          -members=>\@mem);

  #access the family

  foreach my $mem ($family->get_members){
    print $mem->display_id."\t".$mem->desc."\n";
  }

  #select members if members have a Bio::Species Object

  my @mem = $family->get_members(-binomial=>"Homo sapiens");
  @mem = $family->get_members(-ncbi_taxid => 9606);
  @mem = $family->get_members(-common_name=>"Human");
  @mem = $family->get_members(-species=>"sapiens");
  @mem = $family->get_members(-genus=>"Homo");

=head1 DESCRIPTION

This is a simple Family object that may hold any group of object. For more
specific families, one should derive from FamilyI.

=head1 FEEDBACK

Email bioperl-l@bioperl.org for support and feedback.

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::SequenceFamily;

use strict;
use warnings;
use base qw(Bio::Root::Root Bio::Cluster::FamilyI);

=head2 new

 Title   : new
 Usage   : my $family = Bio::Cluster::SequenceFamily->new(
                             -family_id=>"Family_1",
                             -description=>"Family Description Here",
                             -annotation_score=>"100",
                             -members=>\@mem);
 Function: Constructor for SequenceFamily object
 Returns : Bio::Cluster::SequenceFamily object

See L<Bio::Cluster::SequenceFamily>.

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($id,$description,$version,$annot_score,
    $family_score,$members) = $self->_rearrange([qw(FAMILY_ID DESCRIPTION VERSION
                                                    ANNOTATION_SCORE
                                                    FAMILY_SCORE MEMBERS)],@args);
    $self->{'_members'} = [];
    $id             && $self->family_id($id);
    $description    && $self->description($description);
    $version        && $self->version($version);
    $annot_score    && $self->annotation_score($annot_score);
    $family_score   && $self->family_score($family_score);
    $members        && $self->add_members($members);

    return $self;
}

=head2 version

 Title   : version
 Usage   : $family->version("1.0");
 Function: get/set for version
 Returns : a string version of the family generated.

=cut

sub version{
    my ($self,$value) = @_;
    if($value){
        $self->{'_version'} =$value;
    }
    return $self->{'_version'};
}

=head2 annotation_score

 Title   : annotation_score
 Usage   : $family->annotation_score(100);
 Function: get/set for annotation_score which
           represent the confidence in which the
           consensus description has been assigned
           to the family.
 Returns : Bio::SimpleAlign

See L<Bio::SimpleAlign>

=cut

sub annotation_score{
    my ($self,$score) = @_;
    if($score){
        $self->{'_annotation_score'} = $score;
    }
    return $self->{'_annotation_score'};
}

=head2 alignment

 Title   : alignment
 Usage   : $family->alignment($align);
 Function: get/set for an alignment object representing
           the multiple alignment of the members of the family.
 Returns : Bio::SimpleAlign

See L<Bio::SimpleAlign>

=cut

sub alignment {
    my ($self,$align) = @_;
    if($align){
        $self->{'_alignment'} = $align;
    }
    return $self->{'_alignment'};
}

=head2 tree

 Title   : tree
 Usage   : $family->tree($tree);
 Function: get/set for an tree object representing
           the phylogenetic tree of the family.
 Returns : Bio::Tree

See L<Bio::Tree>

=cut

sub tree {
    my ($self,$tree) = @_;
    if($tree) {
        $self->{'_tree'} = $tree;
    }
    return $self->{'_tree'};
}

=head1 L<Bio::Cluster::FamilyI> methods

=cut

=head2 family_score

 Title   : family_score
 Usage   : Bio::Cluster::FamilyI->family_score(95);
 Function: get/set for the score of algorithm used to generate
           the family if present

           This is aliased to cluster_score().

 Returns : the score
 Args    : the score

=cut

sub family_score {
    return shift->cluster_score(@_);
}


=head2 family_id

 Title   : family_id
 Usage   : $family->family_id("Family_1");
 Function: get/set for family id

           This is aliased to display_id().

 Returns : a string specifying identifier of the family

=cut

sub family_id{
    return shift->display_id(@_);
}

=head1 L<Bio::ClusterI> methods

=cut

=head2 display_id

 Title   : display_id
 Usage   :
 Function: Get/set the display name or identifier for the cluster
 Returns : a string
 Args    : optional, on set the display ID ( a string)

=cut

sub display_id{
    my ($self,$id) = @_;
    if($id){
        $self->{'_cluster_id'} = $id;
    }
    return $self->{'_cluster_id'};
}

=head2 description

 Title   : description
 Usage   : $fam->description("POLYUBIQUITIN")
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

=head2 get_members

 Title   : get_members
 Usage   : Valid criteria:
           -common_name
           -binomial
           -ncbi_taxid
           -organelle
           -genus
           $family->get_members(-common_name =>"human");
           $family->get_members(-species     =>"homo sapiens");
           $family->get_members(-ncbi_taxid  => 9606);
           For now, multiple critieria are ORed.

           Will return all members if no criteria are provided.

 Function: get members using methods from L<Bio::Species>
           the phylogenetic tree of the family.
 Returns : an array of objects that are member of this family.

=cut

sub get_members {
    my $self = shift;
    return @{$self->{'_members'}} unless @_;

    ## since the logic behind the checks is OR, we keep the ids in an hash for
    ## performance (skip the test if it's already there) and to avoid repats
    my %match;
    my %filter = @_;
    foreach my $key (keys %filter) {
        (my $method = $key) =~ s/^-//;
        %match = (%match, map { $_ => $_ } grep {
            ! $match{$_} && $_->species &&
            ($_->species->can($method) ||
                $self->throw("$method is an invalid criteria")) &&
            $_->species->$method() eq $filter{$key}
        } @{$self->{'_members'}});
    }
    return map {$match{$_}} keys (%match);
}

=head2 size

 Title   : size
 Usage   : $fam->size();
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
 Usage   : $fam->cluster_score(100);
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


=head1 Implementation specific methods

  These are mostly for adding/removing/changing.

=cut

=head2 add_members

 Title   : add_members
 Usage   : $fam->add_member([$seq1,$seq1]);
 Function: add members to a family
 Returns :
 Args    : the member(s) to add, as an array or arrayref

=cut

sub add_members{
    my ($self,@mems) = @_;

    if (@mems) {
        my $mem = shift(@mems);
        if(ref($mem) eq "ARRAY"){
            push @{$self->{'_members'}},@{$mem};
        } else {
            push @{$self->{'_members'}},$mem;
        }
        push @{$self->{'_members'}}, @mems;
    }
    return 1;
}

=head2 remove_members

 Title   : remove_members
 Usage   : $fam->remove_members();
 Function: remove all members from a family
 Returns : the previous array of members
 Args    : none

=cut

sub remove_members{
    my ($self) =  @_;
    my $mems = $self->{'_members'};
    $self->{'_members'} = [];
    return @$mems;
}

#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*flush_members = \&remove_members;
*add_member = \&add_members;

=head2 members

 Title   : members
 Usage   : $members = $fam->members([$seq1,$seq1]);
 Function: Deprecated. Use add_members() or get_members() instead

=cut

sub members{
    my $self = shift;
    if(@_) {
        # this is in set mode
        $self->warn("setting members() in ".ref($self)." is deprecated.\n".
                    "Use add_members() instead.");
        return $self->add_members(@_);
    } else {
        # get mode
        $self->warn("members() in ".ref($self)." is deprecated.\n".
                    "Use get_members() instead.");
        return $self->get_members();
    }
}

1;
