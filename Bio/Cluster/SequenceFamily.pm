#
# BioPerl module for Bio::Cluster::SequenceFamily
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

  use Bio::Cluster::SequenceFamily

  use Bio::SeqIO;
  use Bio::Cluster::SequenceFamily;

  my $file =  Bio::Root::IO->catfile('t','data','swiss.dat');
  my $seqio= new Bio::SeqIO('-format' => 'swiss',
                          '-file'   => $file);
  my @mem;
  while(my $seq = $seqio->next_seq){
    push @mem, $seq;
  }

  #create the family
  my $family = Bio::Cluster::SequenceFamily->new(-family_id=>"Family_1",
                                                 -description=>"Family Description Here",
                                                 -annotation_score=>"100",
                                               -members=>\@mem);

  #access the family

  foreach my $mem ($family->members){
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


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Shawn Hoon 

Email shawnh@fugu-sg.org 


=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::SequenceFamily;

use strict;
use vars qw(@ISA);


use Bio::Root::Root;
use Bio::Cluster::FamilyI;

@ISA = qw(Bio::Root::Root Bio::Cluster::FamilyI);


=head2 new

 Title   : new
 Usage   : my $family = Bio::Cluster::SequenceFamily->new(-family_id=>"Family_1",
                                       -description=>"Family Description Here",
                                       -annotation_score=>"100",
                                       -members=>\@mem);
 Function: Constructor for SequenceFamily object
 Returns : L<Bio::Cluster::SequenceFamily> object

=cut

sub new {
	my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($id,$description,$version,$annot_score,
  $family_score,$members) = $self->_rearrange([qw(FAMILY_ID DESCRIPTION VERSION 
                                                   ANNOTATION_SCORE 
                                                   FAMILY_SCORE MEMBERS)],@args);
  $id && $self->family_id($id);
  $description && $self->description($description);
  $version && $self->version($version);
  $annot_score && $self->annotation_score($annot_score);
  $family_score && $self->family_score($family_score);
  $members && $self->members($members);

  return $self;

}


######get/sets here##############

=head2 family_id

 Title   : family_id
 Usage   : $family->family_id("Family_1"); 
 Function: get/set for family id 
 Returns : a string specifying identifier of the family 

=cut

sub family_id{
  my ($self,$id) = @_;
  if($id){
    $self->{'_family_id'} = $id;
  }
  return $self->{'_family_id'};
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
 Returns : L<Bio::SimpleAlign> 

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
 Returns : L<Bio::SimpleAlign>

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
 Returns : L<Bio::Tree> 

=cut

sub tree {
  my ($self,$tree) = @_;
  if($tree) {
    $self->{'_tree'} = $tree;
  }
  return $self->{'_tree'};
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
 Function: get members using methods from L<Bio::Species>
           the phylogenetic tree of the family.
 Returns : an array of objects that are member of this family. 

=cut

sub get_members {
    my ($self,@args) = @_;
    my %hash = @args;

    my @return;
    foreach my $mem($self->members){
      foreach my $key (keys %hash){
        my $method = $key;
        $method=~s/-//g;
        if($mem->can('species')){
          my $species = $mem->species;
          $species->can($method) || $self->throw("$method is an invalid criteria");
          if($species->$method eq $hash{$key}){
            push @return, $mem;
          }
        }
      }
    }

    return @return;
}

1;
