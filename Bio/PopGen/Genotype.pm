#
# BioPerl module for Bio::PopGen::Genotype
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

Bio::PopGen::Genotype - An implementation of GenotypeI which is just an allele container

=head1 SYNOPSIS

  use Bio::PopGen::Genotype;
  my $genotype = Bio::PopGen::Genotype->new(-marker_name   => $name,
                                           -individual_id => $indid,
                                           -alleles       => \@alleles);

=head1 DESCRIPTION

This object will contain alleles for a given marker for a given
individual.

The class variable BlankAlleles (accessible through
$Bio::PopGen::Genotype::BlankAlleles = 'somepattern') can be set to a
regexp pattern for identifying blank alleles which should no be
counted (they are effectively missing data).  By default it set to
match white space, '-', 'N' or 'n', and '?' as blank alleles which are
skipped.

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


package Bio::PopGen::Genotype;
use vars qw($BlankAlleles);
use strict;

$BlankAlleles = '[\s\-Nn\?]';


# Object preamble - inherits from Bio::Root::Root



use base qw(Bio::Root::Root Bio::PopGen::GenotypeI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Genotype->new();
 Function: Builds a new Bio::PopGen::Genotype object 
 Returns : an instance of Bio::PopGen::Genotype
 Args    : -marker_name   => string representing name of the marker
           -individual_id => string representing individual id (optional)
           -alleles       => arrayref with each item in the array being an allele

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($marker_name, $marker_type, $ind_id, $alleles) = $self->_rearrange([qw(MARKER_NAME
                                                               MARKER_TYPE
							       INDIVIDUAL_ID
							       ALLELES)],@args);
  defined $marker_name && $self->marker_name($marker_name);
  defined $marker_type && $self->marker_type($marker_type);
  defined $ind_id      && $self->individual_id($ind_id);
  if( defined $alleles ) {
      if( ref($alleles) =~ /array/i ) {
	  $self->add_Allele(@$alleles);
      } else { 
	  $self->warn("Could not initialize with -alleles value, it is not an array ref");
      }
  }
  return $self;
}


=head2 marker_name

 Title   : marker_name
 Usage   : my $name = $genotype->marker_name();
 Function: Get the marker name for a genotype result
 Returns : string
 Args    : [optional] marker name value to store


=cut

sub marker_name{
    my ($self) = shift;
    return $self->{'_marker_name'} = shift if @_;
    return $self->{'_marker_name'};
}

=head2 marker_type

 Title   : marker_type
 Usage   : my $name = $genotype->marker_type();
 Function: Get the marker type for a genotype result
 Returns : M (microsatellite, or other multi-allelic 
           locus) or S (biallelic/SNP locus)
 Args    : [optional] marker type value to store


=cut

sub marker_type{
    my ($self) = shift;
    return $self->{'_marker_type'} = shift if @_;
    return $self->{'_marker_type'};
}


=head2 individual_id

 Title   : individual_id
 Usage   : my $indid = $genotype->individual_id();
 Function: Gets the individual id associated with a genotype
           This is effectively a back reference since we will typically
           associate a genotype with an individual with an 
           individual HAS-A genotype relationship.
 Returns : unique id string for an individual
 Args    : none


=cut

sub individual_id {
    my ($self) = shift;
    return $self->{'_individual_id'} = shift if @_;
    return $self->{'_individual_id'};
}

=head2 get_Alleles

 Title   : get_Alleles
 Usage   : my @alleles = $genotype->get_Alleles();
 Function: Get the alleles for a given marker and individual
 Returns : array of alleles (strings in this implementation)
 Args    : $showblank - boolean flag to indicate return ALL alleles not 
                        skipping the coded EMPTY alleles

 Note    : Uses the class variable $BlankAlleles to test if alleles
           should be skipped or not.

=cut

sub get_Alleles{
    my ($self) = shift;
    
     if( @_ && $_[0] ) {
	return @{$self->{'_alleles'} || []};
    } else {
	if( defined $self->{'_cached_noblank'} ) {
	    return @{$self->{'_cached_noblank'}} 
	}
	# one liners - woo hoo.
	$self->{'_cached_noblank'} = [ grep { ! /^\s*$BlankAlleles\s*$/o } 
				       @{$self->{'_alleles'} || []}];
	return @{$self->{'_cached_noblank'}};
    }
}

=head2 add_Allele

 Title   : add_Allele
 Usage   : $genotype->add_Allele(@alleles);
 Function: Add alleles to the genotype, at this point there is no
           verification to insure that haploid individuals only have 1 
           allele or that diploids only have 2 - we assume that is
           done by the user creating these objects 
 Returns : count of the number of alleles in genotype
 Args    : Array of alleles to store


=cut

sub add_Allele {
    my ($self) = shift;
    $self->{'_cached_noblank'} = undef;    
    push @{$self->{'_alleles'}}, @_;
    return scalar @{$self->{'_alleles'}};
}

=head2 reset_Alleles

 Title   : reset_Alleles
 Usage   : $genotype->reset_Alleles;
 Function: Resets the stored alleles so the list is empty
 Returns : None
 Args    : None


=cut

sub reset_Alleles{
   my ($self,@args) = @_;
   $self->{'_cached_noblank'} = undef;
   $self->{'_alleles'} = [];
   return 0;
}


1;
