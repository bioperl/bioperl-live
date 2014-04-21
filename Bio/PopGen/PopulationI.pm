#
# BioPerl module for Bio::PopGen::PopulationI
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

Bio::PopGen::PopulationI - Interface for Populations

=head1 SYNOPSIS

  # Get Bio::PopGen::PopulationI object somehow, like
  # from Bio::Population::Population

  print "name is ", $population->name(), "\n";
  print "source is ", $population->source(), "\n";
  print "description is ", $population->description(), "\n";

  print "For marker $markername:\n";
  foreach my $genotype ( $population->get_Genotypes(-marker => $markername) ) {
      print "Individual ", $genotype->individual_id, " genotype alleles are ",
      join(',', $genotype->get_Alleles()), "\n";
  }
  # get a marker with allele frequencies calculated from the population
  my $marker = $population->get_Marker($markername); 
  my %af = $marker->get_Allele_Frequencies;
  foreach my $allele ( keys %af ) {
      print "$allele $af{$allele}\n";
  }

=head1 DESCRIPTION

This interface describes the basics of a population.  One can use this
object to get the genotypes of specific individuals, only those
individuals which have a certain marker, or create a marker with
allele frequency information.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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


package Bio::PopGen::PopulationI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

=head2 name

 Title   : name
 Usage   : my $name = $pop->name
 Function: Get the population name
 Returns : string representing population name
 Args    : [optional] string representing population name


=cut

sub name{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}



=head2 description

 Title   : description
 Usage   : my $description = $pop->description
 Function: Get the population description
 Returns : string representing population description
 Args    : [optional] string representing population description


=cut

sub description{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 source

 Title   : source
 Usage   : my $source = $pop->source
 Function: Get the population source
 Returns : string representing population source
 Args    : [optional] string representing population source


=cut

sub source{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 annotation

 Title   : annotation
 Usage   : my $annotation_collection = $pop->annotation;
 Function: Get/set a Bio::AnnotationCollectionI for this population
 Returns : Bio::AnnotationCollectionI object
 Args    : [optional set] Bio::AnnotationCollectionI object


=cut

sub annotation{
   my ($self) = @_;
   $self->throw_not_implemented();
}

=head2 get_Individuals

 Title   : get_Individuals
 Usage   : my @inds = $pop->get_Individuals();
 Function: Return the individuals, alternatively restrict by a criteria
 Returns : Array of L<Bio::PopGen::IndividualI> objects
 Args    : none if want all the individuals OR,
           -unique_id => To get an individual with a specific id
           -marker    => To only get individuals which have a genotype specific
                        for a specific marker name


=cut


sub get_Individuals{
    shift->throw_not_implemented();
}

=head2 get_Genotypes

 Title   : get_Genotypes
 Usage   : my @genotypes = $pop->get_Genotypes(-marker => $name)
 Function: Get the genotypes for all the individuals for a specific
           marker name
 Returns : Array of L<Bio::PopGen::GenotypeI> objects
 Args    : -marker => name of the marker


=cut

sub get_Genotypes{
    shift->throw_not_implemented;
}

=head2 get_Marker

 Title   : get_Marker
 Usage   : my $marker = $population->get_Marker($name)
 Function: Get a Bio::PopGen::Marker object based on this population
 Returns : L<Bio::PopGen::MarkerI> object
 Args    : name of the marker


=cut

sub get_Marker{
    shift->throw_not_implemented();
}

=head2 get_marker_names

 Title   : get_marker_names
 Usage   : my @names = $pop->get_marker_names;
 Function: Get the names of the markers
 Returns : Array of strings
 Args    : none


=cut

sub get_marker_names{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 get_Markers

 Title   : get_Markers
 Usage   : my @markers = $pop->get_Markers();
 Function: Will retrieve a list of instantiated MarkerI objects 
           for a population.  This is a convience method combining
           get_marker_names with get_Marker
 Returns : List of array of Bio::PopGen::MarkerI objects
 Args    : none


=cut

sub get_Markers{
    my ($self) = shift;
    return map { $self->get_Marker($_) } $self->get_marker_names();
}


=head2 get_number_individuals

 Title   : get_number_individuals
 Usage   : my $count = $pop->get_number_individuals;
 Function: Get the count of the number of individuals
 Returns : integer >= 0
 Args    : [optional] marker name, will return a count of the number
           of individuals which have this marker


=cut

sub get_number_individuals{
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;
