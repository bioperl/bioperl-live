
#
# BioPerl module for Bio::Pfam::Annotation::Comment
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Species - Generic species object

=head1 SYNOPSIS

    $species = Bio::Species->new(); # Can also pass classification
                                    # array to new as below
                                    
    $species->classification(qw( sapiens Homo Hominidae
                                 Catarrhini Primates Eutheria
                                 Mammalia Vertebrata Chordata
                                 Metazoa Eukaryota ));
    
    $genus = $species->genus();
    
    $bi = $species->binomial();     # $bi is now "Homo sapiens"
    
    # For storing common name
    $species->common_name("human");

=head1 DESCRIPTION

Provides a very simple object for storing phylogenetic
information.  The classification is stored in an array,
which is a list of nodes in a phylogenetic tree.  Access to
getting and setting species and genus is provided, but not
to any of the other node types (eg: "phlum", "class",
"order", "family").  There's plenty of scope for making the
model more sophisticated, if this is ever needed.

A method is also provided for storing a common name of the
species.

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

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'classification'} = [];
  $self->{'common_name'} = undef;
  if (@args) {
    $self->classification(@args);
  }
  return $make; # success - we hope!
}

=head2 classification

 Title   : classification
 Usage   : $self->classification(@class_array);
           @classification = $self->classification();
 Function: Fills or returns the classifcation list in
           the object.  The array provided must be in
           the order SUBSPECIES, SPECIES, GENUS ---> KINGDOM.
           The first and second element of the array, the subspecies and
           species, must be in lower case, and the rest in title
           case.  Only species must be present.

           Note that the format convention given above has changed after 
           release 0.60. Formerly, SUBSPECIES was not necessary. In order to
           break as few scripts as possible, the method tries to recognize
           whether or not the subspecies is provided, given that the rest
           is given in correct case. This is the reason that the example given
           below is still valid.
 Example : $obj->classification(qw( sapiens Homo Hominidae
           Catarrhini Primates Eutheria Mammalia Vertebrata
           Chordata Metazoa Eukaryota));
 Returns : Classification array
 Args    : Classification array

=cut



sub classification {
    my $self = shift;

    if (@_) {
        my @classification = @_;
        
	# Try to be smart and check whether someone may have omitted the
	# supspecies, but provided everything else in correct case spelling.
	# Likewise, providing the subspecies alone obviously doesn't make
	# sense, so it's probably the species.
	# HL <Hilmar.Lapp@pharma.novartis.com>
	if(scalar(@classification) == 0) {
	    $self->throw("no elements in classification: " .
			 "must at least supply species name");
	} elsif((scalar(@classification) == 1) ||
		(($classification[0] =~ /^[a-z].*/) &&
		 ($classification[1] =~ /^[A-Z][a-z].*/))) {
	    # looks like subspecies has been omitted, so add it
	    unshift(@classification, '');
	}
        # Check the names supplied in the classification string
        {
            # Species should be in lower case
            my $species = $classification[1];
            $self->validate_species_name( $species );

            # All other names must be in title case
            for (my $i = 2; $i < @classification; $i++) {
                $self->validate_name( $classification[$i] );
            }
        }

        # Store classification
        $self->{'classification'} = [ @classification ];
    }
    return @{$self->{'classification'}};
}

=head2 

 Title   : common_name
 Usage   : $self->common_name( $common_name );
           $common_name = $self->common_name();
 Function: Get or set the commonn name of the species
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
        $self->{'classification'}[1] = $species;
    } else {
        return $self->{'classification'}[1];
    }
}

=head2 sub_species

 Title   : sub_species
 Usage   : $obj->sub_species($newval)
 Function: 
 Returns : value of sub_species
 Args    : newvalue (optional)


=cut

sub sub_species{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'classification'}[0] = $value;
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
        $self->{'classification'}[2] = $genus;
    } else {
        return $self->{'classification'}[2];
    }

}

=head2 binomial

 Title   : binomial
 Usage   : $binomial = $self->binomial();
           $binomial = $self->binomial('FULL');
 Function: Returns a string "Genus species", or "Genus species subspecies",
           the first argument is 'FULL'.
 Args    : Optionally the string 'FULL' to get the full name including the
           the subspecies.

=cut


sub binomial {
    my( $self, $full ) = @_;
    
    my( $ssp, $species, $genus ) = $self->classification();
    if(defined($full) && ($full eq 'FULL')) {
	return "$genus $species $ssp";
    } else {
	return "$genus $species";
    }
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


1;

__END__
