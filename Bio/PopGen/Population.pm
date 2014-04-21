#
# BioPerl module for Bio::PopGen::Population
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

Bio::PopGen::Population - A population of individuals

=head1 SYNOPSIS

  use Bio::PopGen::Population;
  use Bio::PopGen::Individual;
  my $population = Bio::PopGen::Population->new();
  my $ind = Bio::PopGen::Individual->new(-unique_id => 'id');
  $population->add_Individual($ind);

  for my $ind ( $population->get_Individuals ) {
    # iterate through the individuals
  }

  for my $name ( $population->get_marker_names ) {
    my $marker = $population->get_Marker($name);
  }

  my $num_inds = $population->get_number_individuals;

  my $homozygote_f   = $population->get_Frequency_Homozygotes;
  my $heterozygote_f = $population->get_Frequency_Heterozygotes;

  # make a population haploid by making fake chromosomes through
  # haplotypes -- ala allele 1 is on chrom 1 and allele 2 is on chrom 2 
  # the number of individuals created will thus be 2 x number in
  # population
  my $happop = $population->haploid_population;


=head1 DESCRIPTION

This is a collection of individuals.  We'll have ways of generating
L<Bio::PopGen::MarkerI> objects out so we can calculate allele_frequencies
for implementing the various statistical tests.

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


package Bio::PopGen::Population;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::PopGen::Marker;
use Bio::PopGen::Genotype;
our $CheckISA = 1;
use base qw(Bio::Root::Root Bio::PopGen::PopulationI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Population->new();
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
      $inds,$checkisa) = $self->_rearrange([qw(NAME 
				     SOURCE 
				     DESCRIPTION
				     INDIVIDUALS
				     CHECKISA)], @args);
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
  $self->{'_checkisa'} = defined $checkisa ? $checkisa : $CheckISA;
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

=head2 annotation

 Title   : annotation
 Usage   : my $annotation_collection = $pop->annotation;
 Function: Get/set a Bio::AnnotationCollectionI for this population
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

=head2 set_Allele_Frequency

 Title   : set_Allele_Frequency
 Usage   : $population->set_Allele_Frequency('marker' => { 'allele1' => 0.1});
 Function: Sets an allele frequency for a Marker for this Population
           This allows the Population to not have individual individual
           genotypes but rather a set of overall allele frequencies
 Returns : Count of the number of markers
 Args    : -name      => (string) marker name
           -allele    => (string) allele name
           -frequency => (double) allele frequency - must be between 0 and 1
           OR
	   -frequencies => { 'marker1' => { 'allele1' => 0.01,
					    'allele2' => 0.99},
			     'marker2' => ...
			    }

=cut

sub set_Allele_Frequency {
   my ($self,@args) = @_;
   my ($name,$allele, $frequency,
       $frequencies) = $self->_rearrange([qw(NAME
					     ALLELE
					     FREQUENCY
					     FREQUENCIES
					     )], @args);
   if( defined $frequencies ) { # this supercedes the res
       if( ref($frequencies) =~ /HASH/i ) {
	   my ($markername,$alleles);
	   while( ($markername,$alleles) = each %$frequencies ) {
	       $self->{'_allele_freqs'}->{$markername} = 
		   Bio::PopGen::Marker->new(-name        => $markername,
					   -allele_freq => $alleles);
	   }
       } else { 
	   $self->throw("Must provide a valid hashref for the -frequencies option");
       }
   } else { 
       unless( defined $self->{'_allele_freqs'}->{$name} ) {
	   $self->{'_allele_freqs'}->{$name} = 
	       Bio::PopGen::Marker->new(-name        => $name);
       }
       $self->{'_allele_freqs'}->{$name}->add_Allele_Frequency($allele,$frequency);
   }
   return scalar keys %{$self->{'_allele_freqs'}};
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
	
	unless( $self->{'_checkisa'} ? $i->isa('Bio::PopGen::IndividualI') : 1  ) {
	    $self->warn("cannot add an individual ($i) which is not a Bio::PopGen::IndividualI");
	    next;
	}
    }
    push @{$self->{'_individuals'}}, @inds;
    $self->{'_cached_markernames'} = undef;
    $self->{'_allele_freqs'} = {};
    return scalar @{$self->{'_individuals'} || []};
}


=head2 remove_Individuals

 Title   : remove_Individuals
 Usage   : $population->remove_Individuals(@ids);
 Function: Remove individual(s) to a population
 Returns : count of the current number in the object 
 Args    : Array of ids

=cut

sub remove_Individuals {
    my ($self,@names) = @_;
    my $i = 0;
    my %namehash; # O(1) lookup will be faster I think
    foreach my $n ( @names ) { $namehash{$n}++ }
    my @tosplice;
    foreach my $ind (  @{$self->{'_individuals'} || []} ) {
	unshift @tosplice, $i if( $namehash{$ind->unique_id} );
	$i++;
    }
    foreach my $index ( @tosplice ) {
	splice(@{$self->{'_individuals'}}, $index,1);
    }
    $self->{'_cached_markernames'} = undef;
    $self->{'_allele_freqs'} = {};
    return scalar @{$self->{'_individuals'} || []};
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
   my @inds = @{$self->{'_individuals'} || []};
   return unless @inds;
   if( @args ) { # save a little time here if @args is empty
       my ($id,$marker) = $self->_rearrange([qw(UNIQUE_ID MARKER)], @args);

       
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
       @{$self->{'_individuals'} || []}
   } 
   $self->warn("You needed to have provided a valid -marker value");
   return ();
}


=head2 get_marker_names

 Title   : get_marker_names
 Usage   : my @names = $pop->get_marker_names;
 Function: Get the names of the markers
 Returns : Array of strings
 Args    : [optional] boolean flag to ignore internal cache status


=cut

sub get_marker_names {
    my ($self,$force) = @_;
    return @{$self->{'_cached_markernames'} || []} 
      if( ! $force && defined $self->{'_cached_markernames'});
    my %unique;
    foreach my $n ( map { $_->get_marker_names } $self->get_Individuals() ) {
	$unique{$n}++;
    }
    my @nms = keys %unique;
    if( $nms[0] =~ /^(Site|Codon)/ ) {
	# sort by site or codon number and do it in 
	# a schwartzian transformation baby!
	@nms = map { $_->[1] } 
 	       sort { $a->[0] <=> $b->[0] }
	       map { [$_ =~ /^(?:Codon|Site)-(\d+)/, $_] } @nms;
    }
    $self->{'_cached_markernames'} = [ @nms ];
    return @{$self->{'_cached_markernames'} || []};
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
   my $marker;
   # setup some caching too
   if( defined $self->{'_allele_freqs'} &&
       defined ($marker = $self->{'_allele_freqs'}->{$markername}) ) {
       # marker is now set to the stored value
   } else { 
       my @genotypes = $self->get_Genotypes(-marker => $markername);
       $marker = Bio::PopGen::Marker->new(-name   => $markername);

       if( ! @genotypes ) {
	   $self->warn("No genotypes for Marker $markername in the population");
       } else { 
	   my %alleles;
	   my $count;
	   for my $al ( map { $_->get_Alleles} @genotypes ) {
	       next if($al eq '?');
	       $count++; 
	       $alleles{$al}++
	   }
	   foreach my $allele ( keys %alleles ) {
	       $marker->add_Allele_Frequency($allele, $alleles{$allele}/$count);
	       $marker->{_marker_coverage} = $count/2;
	   }
       }
       $self->{'_allele_freqs'}->{$markername} = $marker;
   }
   return $marker;
}


=head2 get_number_individuals

 Title   : get_number_individuals
 Usage   : my $count = $pop->get_number_individuals;
 Function: Get the count of the number of individuals
 Returns : integer >= 0
 Args    : none


=cut

sub get_number_individuals{
   my ($self,$markername) = @_;

   if( $self->{'_forced_set_individuals'} ) {
       return $self->{'_forced_set_individuals'};
   }

   unless( defined $markername ) {
       return scalar @{$self->{'_individuals'} || []};
   } else { 
       my $number =0;
       foreach my $individual ( @{$self->{'_individuals'} || []} ) {
	   $number++ if( $individual->has_Marker($markername));
       }
       return $number;
   }
}

=head2 set_number_individuals

 Title   : set_number_individuals
    Usage   : $pop->set_number_individuals($num);
 Function: Fixes the number of individuals, call this with
           0 to unset.
           Only use this if you know what you are doing,
           this is only relavent when you are just adding
           allele frequency data for a population and want to
           calculate something like theta
 Returns : none
 Args    : individual count, calling it with undef or 0
            will reset the value to return a number
            calculated from the number of individuals
            stored for this population.

=cut

sub set_number_individuals{
   my ($self,$indcount) = @_;
   return $self->{'_forced_set_individuals'} = $indcount;
}


=head2 get_Frequency_Homozygotes

 Title   : get_Frequency_Homozygotes
 Usage   : my $freq = $pop->get_Frequency_Homozygotes;
 Function: Calculate the frequency of homozygotes in the population
 Returns : fraction between 0 and 1
 Args    : $markername


=cut

sub get_Frequency_Homozygotes{
   my ($self,$marker,$allelename) = @_;
   my ($homozygote_count) = 0;
   return 0 if ! defined $marker || ! defined $allelename;
   $marker = $marker->name if( defined $marker && 
			       ref($marker) &&
			       ( $self->{'_checkisa'} ? 
				 $marker->isa('Bio::PopGen::MarkerI') : 1));
   my $total = $self->get_number_individuals($marker);
   foreach my $genotype ( $self->get_Genotypes($marker) ) {
       my %alleles = map { $_ => 1} $genotype->get_Alleles();
       # what to do for non-diploid situations?
       if( $alleles{$allelename} ) {
	   $homozygote_count++ if( keys %alleles == 1);
       }
   }
   return $total ? $homozygote_count / $total : 0;
}

=head2 get_Frequency_Heterozygotes

 Title   : get_Frequency_Heterozygotes
 Usage   : my $freq = $pop->get_Frequency_Homozygotes;
 Function: Calculate the frequency of homozygotes in the population
 Returns : fraction between 0 and 1
 Args    : $markername


=cut

sub get_Frequency_Heterozygotes{
   my ($self,$marker,$allelename) = @_;
   my ($heterozygote_count) = 0;
   return 0 if ! defined $marker || ! defined $allelename;
   $marker = $marker->name if( defined $marker && ref($marker) &&
			       ($self->{'_checkisa'} ? 
				$marker->isa('Bio::PopGen::MarkerI') : 1));
   if( ref($marker) ) {
       $self->warn("Passed in a ".ref($marker). " to has_Marker, expecting either a string or a Bio::PopGen::MarkerI");
       return 0;
   }
   my $total = $self->get_number_individuals($marker);

   foreach my $genotype ( $self->get_Genotypes($marker) ) {
       my %alleles = map { $_ => 1} $genotype->get_Alleles();
       # what to do for non-diploid situations?
       if( $alleles{$allelename} ) {
	   $heterozygote_count++ if( keys %alleles == 2);
       }
   }
   return $total ? $heterozygote_count / $total : 0;
}

=head2 haploid_population

 Title   : haploid_population
 Usage   : my $pop = $population->haploid_population;
 Function: Make a new population where all the individuals
           are haploid - effectively an individual out of each
           chromosome an individual has.  
 Returns : L<Bio::PopGen::PopulationI>
 Args    : None


=cut

sub haploid_population{
   my ($self) = @_;
   my @inds;
   my @marker_names = $self->get_marker_names;

   for my $ind ( $self->get_Individuals ) {
       my @chromosomes;
       my $id = $ind->unique_id;
       # separate genotypes into 'chromosomes'
       for my $marker_name( @marker_names ) {
	   my ($genotype) = $ind->get_Genotypes(-marker => $marker_name);
	   my $i =0;
	   for my $allele ( $genotype->get_Alleles ) {
	       push @{$chromosomes[$i]}, 
	       Bio::PopGen::Genotype->new(-marker_name => $marker_name,
					-individual_id => $id.".$i",
					-alleles     => [$allele]);
	       $i++;
	   }
       }
       for my $chrom ( @chromosomes ) {
	   my $copyind = ref($ind)->new(-unique_id => $id.".1",
					-genotypes => $chrom);
	   push @inds, $ind;
       }
   }
   my $population = ref($self)->new(-name        => $self->name,
				    -source      => $self->source,
				    -description => $self->description,
				    -individuals => \@inds);
				    
}

1;
