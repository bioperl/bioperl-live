# $Id $
#
# BioPerl module for Bio::PopGen::IndividualI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IndividualI - An individual who has Genotype or Sequence Results

=head1 SYNOPSIS

  # Get a Bio::PopGen::IndividualI somehow
  # test if it has alleles/genotypes for a given marker
  if( $ind->has_marker($markername) ) {
  }
  # get the unique id
  print $ind->unique_id, "\n";

  # get the number of results (genotypes)
  print $ind->num_results;

=head1 DESCRIPTION

Describe the interface here

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
email or the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::IndividualI;
use strict;


use base qw(Bio::Root::RootI);


=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $individual->unique_id
 Function: Unique Identifier
 Returns : string representing unique identifier
 Args    : string


=cut

sub unique_id{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 num_genotypes

 Title   : num_genotypes
 Usage   : my $count = $person->num_results;
 Function: returns the count of the number of Results for a person
 Returns : integer
 Args    : none

=cut

sub num_genotypes { 
    shift->throw_not_implemented;
}

sub num_of_results{ 
    my $self = shift;
    $self->deprecated("num_of_results is deprecated, use num_genotypes instead");
    $self->num_genotypes;
}

=head2 annotation

 Title   : annotation
 Usage   : my $annotation_collection = $ind->annotation;
 Function: Get/set a Bio::AnnotationCollectionI for this individual
 Returns : Bio::AnnotationCollectionI object
 Args    : [optional set] Bio::AnnotationCollectionI object

=cut

sub annotation{
   my ($self, $arg) = @_;
   $self->throw_not_implemented();
}

=head2 get_Genotypes

 Title   : get_Genotypes
 Usage   : my @genotypes = $ind->get_Genotypes(-marker => $markername);
 Function: Get the genotypes for an individual, based on a criteria
 Returns : Array of genotypes
 Args    : either none (return all genotypes) or 
           -marker => name of marker to return (exact match, case matters)


=cut

sub get_Genotypes{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 has_Marker

 Title   : has_Marker
 Usage   : if( $ind->has_Marker($name) ) {}
 Function: Boolean test to see if an Individual has a genotype 
           for a specific marker
 Returns : Boolean (true or false)
 Args    : String representing a marker name


=cut

sub has_Marker{
   my ($self,$name) = @_;
   $self->throw_not_implemented();
}

=head2 get_marker_names

 Title   : get_marker_names
 Usage   : my @names = $individual->get_marker_names;
 Function: Returns the list of known marker names
 Returns : List of strings
 Args    : none


=cut

sub get_marker_names{
   my ($self) = @_;
   $self->throw_not_implemented();
}


1;
