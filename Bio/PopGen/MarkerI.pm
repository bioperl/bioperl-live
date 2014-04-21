# $Id $
#
# BioPerl module for Bio::PopGen::MarkerI
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

Bio::PopGen::MarkerI - A Population Genetic conceptual marker

=head1 SYNOPSIS

# Get a Bio::PopGen::MarkerI somehow - like using a Bio::PopGen::Marker

  my $name = $marker->name();            # marker name
  my $description = $marker->description(); # description
  my $type = $marker->type();            # coded type of the marker
  my $unique_id = $marker->unique_id;    # optional unique ID

  my @alleles = $marker->get_Alleles();  # the known alleles
  my %allele_freqs = $marker->get_Allele_Frequencies(); # keys are marker names
                                         # vals are frequencies
                                         # may change to handle multiple populations


=head1 DESCRIPTION

This is the basic interface for Markers which one can associate
alleles with for calculating Theta and Pi.

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

=head1 CONTRIBUTORS

Matthew Hahn, matthew.hahn-at-duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::MarkerI;
use strict;


use base qw(Bio::Root::RootI Bio::AnnotatableI);


=head2 name

 Title   : name
 Usage   : my $name = $marker->name();
 Function: Get the name of the marker
 Returns : string representing the name of the marker
 Args    : 


=cut

sub name{
   $_[0]->throw_not_implemented();
}


=head2 description

 Title   : description
 Usage   : my $desc = $marker->description
 Function: Get the marker description free text
 Returns : string
 Args    : [optional] string


=cut

sub description{
    $_[0]->throw_not_implemented();
}

=head2 type

 Title   : type
 Usage   : my $type = $marker->type;
 Function: Get coded string for marker type
 Returns : string
 Args    : [optional] string


=cut

sub type{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $marker->unique_id;
 Function: Get the unique marker ID
 Returns : unique ID string
 Args    : [optional ] string


=cut

sub unique_id{
   my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($seq_obj)
 Function: retrieve the attached annotation object
 Returns : Bio::AnnotationCollectionI or none;

See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information. This method comes through extension from
L<Bio::AnnotatableI>.


=cut


sub annotation{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 get_Alleles

 Title   : get_Alleles
 Usage   : my @alleles = $marker->get_Alleles();
 Function: Get the available marker alleles if they are known and stored
 Returns : Array of strings
 Args    : none


=cut

sub get_Alleles{
   my ($self) = @_;
   $self->throw_not_implemented();
}


=head2 get_Allele_Frequencies

 Title   : get_Allele_Frequencies
 Usage   : my %allele_freqs = $marker->get_Allele_Frequencies;
 Function: Get the alleles and their frequency (set relative to
           a given population - you may want to create different
           markers with the same name for different populations
           with this current implementation
 Returns : Associative array (hash) where keys are the names of the alleles
 Args    : none


=cut

sub get_Allele_Frequencies{
   my ($self) = @_;
    $self->throw_not_implemented();
}

1;
