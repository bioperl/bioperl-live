# $Id$
#
# BioPerl module for Bio::Species
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Species - Generic species object

=head1 SYNOPSIS

    $species = Bio::Species->new(-classification => [@classification]);
                                    # Can also pass classification
                                    # array to new as below

    $species->classification(qw( sapiens Homo Hominidae
                                 Catarrhini Primates Eutheria
                                 Mammalia Vertebrata Chordata
                                 Metazoa Eukaryota ));

    $genus = $species->genus();

    $bi = $species->binomial();     # $bi is now "Homo sapiens"

    # For storing common name
    $species->common_name("human");

    # For storing subspecies
    $species->sub_species("accountant");

=head1 DESCRIPTION

Provides a very simple object for storing phylogenetic
information.  The classification is stored in an array,
which is a list of nodes in a phylogenetic tree.  Access to
getting and setting species and genus is provided, but not
to any of the other node types (eg: "phylum", "class",
"order", "family").  There's plenty of scope for making the
model more sophisticated, if this is ever needed.

A methods are also provided for storing common
names, and subspecies.

=head1 CONTACT

James Gilbert email B<jgrg@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


#' Let the code begin...


package Bio::Species;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Root;


@ISA = qw(Bio::Root::Root);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'classification'} = [];
  $self->{'common_name'} = undef;
  my ($classification) = $self->_rearrange([qw(CLASSIFICATION)], @args);
  if( defined $classification &&
      (ref($classification) eq "ARRAY") ) {
      $self->classification(@$classification);
  }
  return $self;
}

=head2 classification

 Title   : classification
 Usage   : $self->classification(@class_array);
           @classification = $self->classification();
 Function: Fills or returns the classification list in
           the object.  The array provided must be in
           the order SPECIES, GENUS ---> KINGDOM.
           Checks are made that species is in lower case,
           and all other elements are in title case.
 Example : $obj->classification(qw( sapiens Homo Hominidae
           Catarrhini Primates Eutheria Mammalia Vertebrata
           Chordata Metazoa Eukaryota));
 Returns : Classification array
 Args    : Classification array 
                 OR
           A reference to the classification array. In the latter case
           if there is a second argument and it evaluates to true,
           names will not be validated.


=cut


sub classification {
    my ($self,@args) = @_;

    if (@args) {

	my ($classif,$force);
	if(ref($args[0])) {
	    $classif = shift(@args);
	    $force = shift(@args);
	} else {
	    $classif = \@args;
	}
	
        # Check the names supplied in the classification string
	# Species should be in lower case
	if(! $force) {
	    $self->validate_species_name($classif->[0]);
	    # All other names must be in title case
	    foreach  (@$classif) {
		$self->validate_name( $_ );
	    }
	}
        # Store classification
        $self->{'classification'} = $classif;
    }
    return @{$self->{'classification'}};
}

=head2 common_name

 Title   : common_name
 Usage   : $self->common_name( $common_name );
           $common_name = $self->common_name();
 Function: Get or set the common name of the species
 Example : $self->common_name('human')
 Returns : The common name in a string
 Args    : String, which is the common name (optional)

=cut

sub common_name{
    my $self = shift;

    return $self->{'common_name'} = shift if @_;
    return $self->{'common_name'};
}

=head2 variant

 Title   : variant
 Usage   : $obj->variant($newval)
 Function: Get/set variant information for this species object (strain,
           isolate, etc).
 Example : 
 Returns : value of variant (a scalar)
 Args    : new value (a scalar or undef, optional)


=cut

sub variant{
    my $self = shift;

    return $self->{'variant'} = shift if @_;
    return $self->{'variant'};
}

=head2 organelle

 Title   : organelle
 Usage   : $self->organelle( $organelle );
           $organelle = $self->organelle();
 Function: Get or set the organelle name
 Example : $self->organelle('Chloroplast')
 Returns : The organelle name in a string
 Args    : String, which is the organelle name

=cut

sub organelle {
    my($self, $name) = @_;

    if ($name) {
        $self->{'organelle'} = $name;
    } else {
        return $self->{'organelle'}
    }
}

=head2 species

 Title   : species
 Usage   : $self->species( $species );
           $species = $self->species();
 Function: Get or set the scientific species name.  The species
           name must be in lower case.
 Example : $self->species( 'sapiens' );
 Returns : Scientific species name as string
 Args    : Scientific species name as string

=cut


sub species {
    my($self, $species) = @_;

    if ($species) {
        $self->validate_species_name( $species );
        $self->{'classification'}[0] = $species;
    }
    return $self->{'classification'}[0];
}

=head2 genus

 Title   : genus
 Usage   : $self->genus( $genus );
           $genus = $self->genus();
 Function: Get or set the scientific genus name.  The genus
           must be in title case.
 Example : $self->genus( 'Homo' );
 Returns : Scientific genus name as string
 Args    : Scientific genus name as string

=cut


sub genus {
    my($self, $genus) = @_;

    if ($genus) {
        $self->validate_name( $genus );
        $self->{'classification'}[1] = $genus;
    }
    return $self->{'classification'}[1];
}

=head2 sub_species

 Title   : sub_species
 Usage   : $obj->sub_species($newval)
 Function:
 Returns : value of sub_species
 Args    : newvalue (optional)


=cut

sub sub_species {
    my( $self, $sub ) = @_;

    if ($sub) {
        $self->{'_sub_species'} = $sub;
    }
    return $self->{'_sub_species'};
}

=head2 binomial

 Title   : binomial
 Usage   : $binomial = $self->binomial();
           $binomial = $self->binomial('FULL');
 Function: Returns a string "Genus species", or "Genus species subspecies",
           the first argument is 'FULL' (and the species has a subspecies).
 Args    : Optionally the string 'FULL' to get the full name including
           the subspecies.

=cut


sub binomial {
    my( $self, $full ) = @_;

    my( $species, $genus ) = $self->classification();
    unless( defined $species) {
	$species = 'sp.';
	$self->warn("classification was not set");
    }
    $genus = ''   unless( defined $genus);
    my $bi = "$genus $species";
    if (defined($full) && ((uc $full) eq 'FULL')) {
	my $ssp = $self->sub_species;
        $bi .= " $ssp" if $ssp;
    }
    return $bi;
}

sub validate_species_name {
    my( $self, $string ) = @_;

    return 1 if $string eq "sp.";
    return 1 if $string =~ /^[a-z][\w\s]+$/i;
    $self->throw("Invalid species name '$string'");
}

sub validate_name {
    return 1; # checking is disabled as there is really not much we can
              # enforce HL 2002/10/03
#     my( $self, $string ) = @_;

#     return 1 if $string =~ /^[\w\s\-\,\.]+$/ or
#         $self->throw("Invalid name '$string'");
}

=head2 ncbi_taxid

 Title   : ncbi_taxid
 Usage   : $obj->ncbi_taxid($newval)
 Function: Get/set the NCBI Taxon ID
 Returns : the NCBI Taxon ID as a string
 Args    : newvalue to set or undef to unset (optional)


=cut

sub ncbi_taxid {
    my $self = shift;

    return $self->{'_ncbi_taxid'} = shift if @_;
    return $self->{'_ncbi_taxid'};
}

1;

__END__
