#
# BioPerl module for Bio::PopGen::Simulation::GeneticDrift
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::Simulation::GeneticDrift - A simple genetic drift simulation

=head1 SYNOPSIS

  use Bio::PopGen::Simulation::GeneticDrift;
  my $sim = Bio::PopGen::Simulation::GeneticDrift->new(-popsize => 40,
						      -alleles => {A => 0.2,
							           B => 0.8});
  for(my $i =0 ;$i < 10; $i++ ) {
    my %f = $sim->next_generation; # get the freqs for each generation
  }

  for(my $i =0 ;$i < 10; $i++ ) {
    # get the allele freqs as part of a Bio::PopGen::Population object
    my $pop = $sim->next_generation('population'); 
  }

=head1 DESCRIPTION

A very simple 1 locus multi-allele random drift module, start with an
initial set of allele frequency and simulate what happens over time.

This isn't really useful for anything in particular yet but will be
built upon.

See Gillespie JH. (1998) "Population Genetics: a Concise guide." The Johns
              Hopkins University Press, Baltimore, USA.  pp.19-47.

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

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Simulation::GeneticDrift;
use strict;

use Bio::PopGen::Population;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Simulation::GeneticDrift->new();
 Function: Builds a new Bio::PopGen::Simulation::GeneticDrift object 
 Returns : an instance of Bio::PopGen::Simulation::GeneticDrift
 Args    : -popsize => starting N
           -haploid => boolean if we should simulate haploids 
           -alleles => arrayref of the allele names
           OR
           -population => L<Bio::PopGen::PopulationI> object to initialize 
                          from some previously defined Population object
                          (or result from a previous simulation)

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($population,
      $popsize, $haploid, $alleles) = $self->_rearrange([qw(POPULATION
							    POPSIZE
							    HAPLOID
							    ALLELES)],@args);
  if( defined $population && ref($population) &&
      $population->isa('Bio::PopGen::PopulationI') ) {
      $self->population_size($population->get_number_individuals || $popsize);
      my %f = $population->get_Allele_Frequencies;
      while( my ($allele,$freq) = each %f ) {
	  $self->add_Allele_Frequency($allele,$freq);
      }
  } else { 
      $self->population_size($popsize);  
  
      if( ! defined $alleles || ref($alleles) !~ /HASH/i  ) {
	  $self->throw("Must provide a valid set of initial allele frequencies to $class as an hashref");
      } 
      while( my ($allele,$freq) = each %$alleles ) {
	  $self->add_Allele_Frequency($allele,$freq);
      }
  }
  unless( $self->validate_Frequencies ) {
      $self->throw("You specified allele frequencies which summed to more than 1");
  }
  return $self;
}


=head2 next_generation

 Title   : next_generation
 Usage   : my %generation = $sim->next_generation
 Function: Get the next generation of allele frequencies based on the current
           generation
 Returns : Hash of allele frequencies
 Args    : 'allelefreqs' or 'population' to get back a hash of allele 
                 frequencies (default) OR a L<Bio::PopGen::Population> object


=cut

sub next_generation{
   my ($self,$rettype) = @_;
   my %initial = $self->get_Allele_Frequencies;
   my $popsize = $self->population_size || 
       $self->throw("Need to have set a valid population size when running the simulation");
   # we're going to construct a mapping of the rational space from 0->1 
   # which will map to a particular allele and be proportional to it
   # frequency
   my ($last,@mapping) = (0);

   # we'll make ranges that cover from >= left and < right in terms of the
   # order doesn't matter - 'distance' does
   # range that we're going to try and match
   # since rand() goes from 0 up to 1 (not including 1)
   foreach my $a ( keys %initial ) {
       push @mapping, [$last,$initial{$a}+$last,$a];
       $last += $initial{$a};
   }

   my %f;
   for( my $i =0; $i < $popsize; $i++ ) {
       my $rand = rand(1);
       foreach my $val ( @mapping ) {
	   if( $rand >= $val->[0] && $rand < $val->[1] ) {
	       $f{$val->[2]}++;
	       last;
	 }
       }
   }
   foreach my $f ( values %f ) {
       $f /= $popsize;
   }
   %{$self->{'_allele_freqs'}} = %f;
   
   if( defined $rettype && 
       $rettype =~ /population/i) {
       return Bio::PopGen::Poulation->new(-frequencies => \%f);
   } else { 
       return %f;
   }

}

=head2 population_size

 Title   : population_size
 Usage   : $obj->population_size($newval)
 Function: 
 Example : 
 Returns : value of population_size (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub population_size{
    my $self = shift;

    return $self->{'_population_size'} = shift if @_;
    return $self->{'_population_size'};
}

=head2 set_Frequencies_Equivalent

 Title   : set_Frequencies_Equivalent
 Usage   : $sim->set_Frequencies_Equivalent
 Function: Reset the allele frequencies so they are all even
 Returns : none
 Args    : none


=cut

sub set_Frequencies_Equivalent{
   my ($self) = @_;
   my @alleles = keys %{$self->{'_allele_freqs'}};
   my $eqfreq  = 1 / scalar @alleles;
   for ( @alleles ) { $self->{'_allele_freqs'}->{$_} = $eqfreq }
   return;
}


=head2 get_Allele_Frequencies

 Title   : get_Allele_Frequencies
 Usage   : my %allele_freqs = $marker->get_Allele_Frequencies;
 Function: Get the alleles and their frequency (set relative to
           a given population - you may want to create different
           markers with the same name for different populations
           with this current implementation
 Returns : Associative array where keys are the names of the alleles
 Args    : none


=cut

sub get_Allele_Frequencies{
   return %{$_[0]->{'_allele_freqs'}};
}

=head2 add_Allele_Frequency

 Title   : add_Allele_Frequency
 Usage   : $marker->add_Allele_Frequency($allele,$freq)
 Function: Adds an allele frequency
 Returns : None
 Args    : $allele - allele name
           $freq   - frequency value


=cut

sub add_Allele_Frequency{
   my ($self,$allele,$freq) = @_;
   $self->{'_allele_freqs'}->{$allele} = $freq;
}

=head2 reset_alleles

 Title   : reset_alleles
 Usage   : $marker->reset_alleles();
 Function: Reset the alleles for a marker
 Returns : None
 Args    : None


=cut

sub reset_alleles{
   my ($self) = @_;
   $self->{'_allele_freqs'} = {};
}

=head2 validate_Frequencies

 Title   : validate_Frequencies
 Usage   : if( $sim->validate_Frequencies) {}
 Function: Sanity checker that allele frequencies sum to 1 or less
 Returns : boolean
 Args    : -strict => 1 boolean if you want to insure that sum of freqs is 1


=cut

sub validate_Frequencies{
   my ($self,@args) = @_;
   my ($strict) = $self->_rearrange([qw(STRICT)], @args);
   my $sum = 0;
   my %freq = $self->get_Allele_Frequencies;
   foreach my $f ( values %freq ) { 
       $sum += $f;
   }
   return ($strict) ? $sum == 1 : $sum <= 1;
}


1;
