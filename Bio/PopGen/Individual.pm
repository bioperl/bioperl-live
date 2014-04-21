#
# BioPerl module for Bio::PopGen::Individual
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

Bio::PopGen::Individual - An implementation of an Individual who has
Genotype or Sequence Results

=head1 SYNOPSIS

  use Bio::PopGen::Individual;

  my $ind = Bio::PopGen::Individual->new(-unique_id => $id,
                                        -genotypes => \@genotypes);

=head1 DESCRIPTION

This object is a container for genotypes.

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
the web:

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


package Bio::PopGen::Individual;
use vars qw($UIDCOUNTER);
use strict;
BEGIN { $UIDCOUNTER = 1 }

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::PopGen::IndividualI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Individual->new();
 Function: Builds a new Bio::PopGen::Individual object 
 Returns : an instance of Bio::PopGen::Individual
 Args    : -unique_id => $id,
           -genotypes => \@genotypes


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->{'_genotypes'} = {};
  my ($uid,$genotypes) = $self->_rearrange([qw(UNIQUE_ID
					       GENOTYPES)],@args);
  unless( defined $uid ) {
      $uid = $UIDCOUNTER++;
  } 
  $self->unique_id($uid);
  if( defined $genotypes ) {
      if( ref($genotypes) =~ /array/i ) {
	  $self->add_Genotype(@$genotypes);
      } else { 
	  $self->warn("Must provide a valid array reference to set the genotypes value in the contructor");
      }
  }
  return $self;
}

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $individual->unique_id
 Function: Unique Identifier
 Returns : string representing unique identifier
 Args    : string


=cut

sub unique_id{
   my ($self) = shift;
   return $self->{'_unique_id'} = shift if @_;
   return $self->{'_unique_id'};
}

=head2 num_of_results

 Title   : num_of_results
 Usage   : my $count = $person->num_results;
 Function: returns the count of the number of Results for a person
 Returns : integer
 Args    : none

=cut

sub num_of_results {
    return scalar keys %{shift->{'_genotypes'}};
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
   return $self->{_annotation} unless $arg;
   $self->throw("Bio::AnnotationCollectionI required for argument") unless
       ref($arg) && $arg->isa('Bio::AnnotationCollectionI');
   return $self->{_annotation} = $arg;
}

=head2 add_Genotype

 Title   : add_Genotype
 Usage   : $individual->add_Genotype
 Function: add a genotype value
 Returns : count of the number of genotypes associated with this individual
 Args    : @genotypes - L<Bio::PopGen::GenotypeI> object(s) containing 
                        alleles plus a marker name

=cut

sub add_Genotype {
   my ($self,@genotypes) = @_;
   
   foreach my $g ( @genotypes ) {
       if( !ref($g) || ! $g->isa('Bio::PopGen::GenotypeI') ) {
	   $self->warn("cannot add $g as a genotype skipping");
	   next;
       }
       my $mname = $g->marker_name;
       if( ! defined $mname || ! length($mname) ) { 
         # can't just say ! name b/c '0' wouldn't be valid 
	   $self->warn("cannot add genotype because marker name is not defined or is an empty string");
	   next;
       }
       if( $self->verbose > 0 && 
	   defined $self->{'_genotypes'}->{$mname} ) {
	   # a warning when we have verbosity cranked up 
	   $self->debug("Overwriting the previous value for $mname for this individual");
       }
       # this will force Genotype individual_id to be set to 
       # the Individual it has been added for
       $g->individual_id($self->unique_id);
       $self->{'_genotypes'}->{$mname} = $g;
   }
   return scalar keys %{$self->{'_genotypes'}};
}

=head2 reset_Genotypes

 Title   : reset_Genotypes
 Usage   : $individual->reset_Genotypes;
 Function: Reset the genotypes stored for this individual
 Returns : none
 Args    : none


=cut

sub reset_Genotypes{
    shift->{'_genotypes'} = {};
}

=head2 remove_Genotype

 Title   : remove_Genotype
 Usage   : $individual->remove_Genotype(@names)
 Function: Removes the genotypes for the requested markers
 Returns : none
 Args    : Names of markers 


=cut

sub remove_Genotype{
   my ($self,@mkrs) = @_;
   foreach my $m ( @mkrs ) {
       delete($self->{'_genotypes'}->{$m});
   }
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
   my ($self,@args) = @_;
   if( @args ) {
       unshift @args, '-marker' if( @args == 1 );  # deal with single args
       
       my ($name) = $self->_rearrange([qw(MARKER)], @args);
       if( ! defined($name) ) {
	   $self->warn("Only know how to process the -marker field currently");
	   return();
       }
       my $v = $self->{'_genotypes'}->{$name};
       return $v;
   }
   return values %{$self->{'_genotypes'} || {}};
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
   return 0 if ! defined $name;

   $name = $name->name if ref($name) && $name->isa('Bio::PopGen::MarkerI');
   if( ref($name) ) { 
       $self->warn("Passed in a ".ref($name). " to has_Marker, expecting either a string or a Bio::PopGen::MarkerI");
       return 0;
   }
   return defined $self->{'_genotypes'}->{$name};
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
   return keys %{$self->{'_genotypes'}};
}


1;
