# $Id$
#
# BioPerl module for Bio::PopGen::Population
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::Population - A population of individuals

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This is a collection of individuals.  We'll have ways of generating
Bio::PopGen::Marker objects out so we can calculate allele_frequencies
for implementing the various statistical tests.

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

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Matthew Hahn <matthew.hahn-at-duke.edu>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Population;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::PopGen::PopulationI;
use Bio::PopGen::Marker;

@ISA = qw(Bio::Root::Root Bio::PopGen::PopulationI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PopGen::Population();
 Function: Builds a new Bio::PopGen::Population object 
 Returns : an instance of Bio::PopGen::Population
 Args    : -individuals => array ref of individuals (optional)
           -name        => population name (optional)
           -source      => a source tag (optional)
           -description => a short description string of the population (optional)

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->{'_individuals'} = [];
  my ($name,$source,$description,
      $inds) = $self->_rearrange([qw(NAME 
				     SOURCE 
				     DESCRIPTION
				     INDIVIDUALS)], @args);
  if( defined $inds ) {
      if( ref($inds) !~ /ARRAY/i ) {
	  $self->warn("Need to provide a value array ref for the -individuals initialization flag");
      } else { 
	  $self->add_Individual(@$inds);
      }
  }

  defined $name   && $self->name($name);
  defined $source && $self->source($source);
  defined $description && $self->description($description);

  return $self;
}


=head2 name

 Title   : name
 Usage   : my $name = $pop->name
 Function: Get the population name
 Returns : string representing population name
 Args    : [optional] string representing population name


=cut

sub name{
   my $self = shift;
   return $self->{'_name'} = shift if @_;
   return $self->{'_name'};
}

=head2 description

 Title   : description
 Usage   : my $description = $pop->description
 Function: Get the population description
 Returns : string representing population description
 Args    : [optional] string representing population description


=cut

sub description{
   my $self = shift;
   return $self->{'_description'} = shift if @_;
   return $self->{'_description'};
}

=head2 source

 Title   : source
 Usage   : my $source = $pop->source
 Function: Get the population source
 Returns : string representing population source
 Args    : [optional] string representing population source


=cut

sub source{
   my $self = shift;
   return $self->{'_source'} = shift if @_;
   return $self->{'_source'};
}

=head2 add_Individual

 Title   : add_Individual
 Usage   : $population->add_Individual(@individuals);
 Function: Add individuals to a population
 Returns : count of the current number in the object 
 Args    : Array of Individuals


=cut

sub add_Individual{
    my ($self,@inds) = @_;
    foreach my $i ( @inds ) {
	next if ! defined $i;
	unless(  $i->isa('Bio::PopGen::IndividualI') ) {
	    $self->warn("cannot add an individual ($i) which is not a Bio::PopGen::IndividualI");
	    next;
	}
	push @{$self->{'_individuals'}}, $i;
    }
    return scalar @{$self->{'_individuals'}};
}

=head2 get_Individuals

 Title   : get_Individuals
 Usage   : my @inds = $pop->get_Individuals();
 Function: Return the individuals, alternatively restrict by a criteria
 Returns : Array of Bio::PopGen::IndividualI objects
 Args    : none if want all the individuals OR,
           -unique_id => To get an individual with a specific id
           -marker    => To only get individuals which have a genotype specific
                        for a specific marker name


=cut

sub get_Individuals{
   my ($self,@args) = @_;
   my @inds = @{$self->{'_individuals'}};
   return unless @inds;
   if( @args ) { # save a little time here if @args is empty
       my ($id,$marker) = $self->_rearrange([qw(UNIQUE_ID
						MARKER)], @args);

       if( defined $id ) { 
	   @inds = grep { $_->unique_id eq $id } @inds;
       } elsif (defined $marker) { 
	   @inds = grep { $_->has_Marker($marker) } @inds;
       }
   }
   return @inds;
}

=head2 get_Genotypes

 Title   : get_Genotypes
 Usage   : my @genotypes = $pop->get_Genotypes(-marker => $name)
 Function: Get the genotypes for all the individuals for a specific
           marker name
 Returns : Array of Bio::PopGen::GenotypeI objects
 Args    : -marker => name of the marker


=cut

sub get_Genotypes{
   my ($self,@args) = @_;
   my ($name) = $self->_rearrange([qw(MARKER)],@args);
   if( defined $name ) {
       return grep { defined $_ } map { $_->get_Genotypes(-marker => $name) } 
       @{$self->{'_individuals'}}
   } 
   $self->warn("You needed to have provided a valid -marker value");
   return ();
}


=head2 get_Marker

 Title   : get_Marker
 Usage   : my $marker = $population->get_Marker($name)
 Function: Get a Bio::PopGen::Marker object based on this population
 Returns : Bio::PopGen::MarkerI object
 Args    : name of the marker


=cut

sub get_Marker{
   my ($self,$markername) = @_;
   my @genotypes = $self->get_Genotypes(-marker => $markername);
   my $marker = new Bio::PopGen::Marker(-name => $markername);
   if( ! @genotypes ) {
       $self->warn("No genotypes for this Marker in the population");
   } else { 
       my %alleles;
       my $count;
       map { $count++; $alleles{$_}++ } map { $_->get_Alleles } @genotypes;
       foreach my $allele ( keys %alleles ) {
	   $marker->add_Allele_Frequency($allele, $alleles{$allele}/$count);
       }
   }
   return $marker;
}




1;
