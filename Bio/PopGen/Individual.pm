# $Id$
#
# BioPerl module for Bio::PopGen::Individual
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::Individual - An implementation of an Individual who has Genotype or Sequence Results

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This object is a container for genotypes

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
the web:

  http://bugzilla.bioperl.org/

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
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::PopGen::IndividualI;

@ISA = qw(Bio::Root::Root Bio::PopGen::IndividualI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PopGen::Individual();
 Function: Builds a new Bio::PopGen::Individual object 
 Returns : an instance of Bio::PopGen::Individual
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->{'_genotypes'} = {};
  my ($uid,$genotypes) = $self->_rearrange([qw(UNIQUE_ID
					       GENOTYPES)],@args);
  defined $uid && $self->unique_id($uid);
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


=head2 add_Genotype

 Title   : add_Genotype
 Usage   : $individual->add_Genotype
 Function: add a genotype value
 Returns : count of the number of genotypes associated with this individual
 Args    : $genotype - Bio::PopGen::GenotypeI object containing the alleles for
                       a marker


=cut

sub add_Genotype {
   my ($self,@genotypes) = @_;
   foreach my $g ( @genotypes ) {
       if( ! defined $g || ! $g->isa('Bio::PopGen::GenotypeI') ) {
	   $self->warn("cannot add genotype, it is not a Bio::PopGen::GenotypeI object");
	   next;
       } elsif( ! length($g->marker_name) ) {
	   $self->warn("cannot add genotype, it must have a valid marker_name associated with it");
       }
       $self->{'_genotypes'}->{$g->marker_name} = $g;
   }
   return scalar keys %{$self->{'_genotypes'}};
}

=head2 reset_Genotypes

 Title   : reset_Genotypes
 Usage   : $genoetype->reset_Genotypes;
 Function: Reset the genotypes stored for this individual
 Returns : none
 Args    : none


=cut

sub reset_Genotypes{
   my ($self,@args) = @_;
   $self->{'_genotypes'} = {};
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
       if( ! $name ) {
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
   return defined $self->{'_genotypes'}->{$name};
}

=head2 get_Marker_Names

 Title   : get_Marker_Names
 Usage   : my @names = $individual->get_Marker_Names;
 Function: Returns the list of known marker names
 Returns : List of strings
 Args    : none


=cut

sub get_Marker_Names{
   my ($self) = @_;
   return keys %{$self->{'_genotypes'}};
}


1;
