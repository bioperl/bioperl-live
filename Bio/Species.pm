
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
      (ref($classification) =~ /array/i) ) {
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

=cut



sub classification {
    my ($self,@args) = @_;

    if (@args) {

        # Check the names supplied in the classification string
        {
            # Species should be in lower case
            my $species = $args[0];
            $self->validate_species_name( $species );

            # All other names must be in title case
            for (my $i = 1; $i < @args; $i++) {
                $self->validate_name( $args[$i] );
            }
        }
        # Store classification
        $self->{'classification'} = [ @args ];
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
 Args    : String, which is the common name

=cut

sub common_name {
    my($self, $name) = @_;

    if ($name) {
        $self->{'common_name'} = $name;
    } else {
        return $self->{'common_name'}
    }
}
=head2

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
 Args    : Optionally the string 'FULL' to get the full name including the
           the subspecies.

=cut


sub binomial {
    my( $self, $full ) = @_;

    my( $species, $genus ) = $self->classification();
    unless( defined $species) {
	$species = '';
	$self->warn("classification was not set");
    }
    $genus = ''   unless( defined $genus);
    my $bi = "$genus $species";
    if (defined($full) && ($full eq 'FULL')) {
	my $ssp = $self->sub_species;
        $bi .= " $ssp" if $ssp;
    }
    return $bi;
}

sub validate_species_name {
    my( $self, $string ) = @_;

    $string =~ /^[\S\d\.]+$||""/ or
        $self->throw("Invalid species name '$string'");
}

sub validate_name {
    my( $self, $string ) = @_;

    return $string =~ /^[A-Z][a-z]+$/ or
        $self->throw("Invalid name '$string' (Wrong case?)");
}

=head2 ncbi_taxid

 Title   : ncbi_taxid
 Usage   : $obj->ncbi_taxid($newval)
 Function:
 Returns : value of ncbi_taxid as string
 Args    : newvalue (optional)


=cut

sub ncbi_taxid {
    my( $self, $sub ) = @_;

    if ($sub) {
        $self->{'_ncbi_taxid'} = $sub;
    }
    return $self->{'_ncbi_taxid'};
}

1;

__END__
