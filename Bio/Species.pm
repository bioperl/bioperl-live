# $Id$
#
# BioPerl module for Bio::Species
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
# Reimplemented by Sendu Bala <bix@sendu.me.uk>
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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

James Gilbert email B<jgrg@sanger.ac.uk>

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' Let the code begin...

package Bio::Species;
use strict;

use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Scalar::Util qw(weaken isweak);
use base qw(Bio::Taxon);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Species->new(-classification => \@class)
 Function: Build a new Species object
 Returns : Bio::Species object
 Args    : -ncbi_taxid     => NCBI taxonomic ID (optional)
           -classification => arrayref of classification

=cut

sub new {
    my($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    my ($org, $sp, $var, $classification) =
        $self->_rearrange([qw(ORGANELLE
					       SUB_SPECIES
					       VARIANT
					       CLASSIFICATION)], @args);
    
    if (defined $classification && ref($classification) eq "ARRAY" && @{$classification}) {
        $self->classification(@$classification);
    }
    else {
        # store a tree on ourselves so we can use Tree methods
        $self->{tree} = new Bio::Tree::Tree();
        
        # some things want to freeze/thaw Bio::Species objects, but
        # _root_cleanup_methods contains a CODE ref, delete it.
        # delete $self->{tree}->{_root_cleanup_methods};
    }
    
    defined $org && $self->organelle($org);
    defined $sp  && $self->sub_species($sp); 
    defined $var && $self->variant($var);
    
    return $self;
}

=head2 classification

 Title   : classification
 Usage   : $self->classification(@class_array);
           @classification = $self->classification();
 Function: Get/set the lineage of this species. The array provided must be in
           the order ... ---> SPECIES, GENUS ---> KINGDOM ---> etc.
 Example : $obj->classification(qw( 'Homo sapiens' Homo Hominidae
           Catarrhini Primates Eutheria Mammalia Vertebrata
           Chordata Metazoa Eukaryota));
 Returns : Classification array
 Args    : Classification array 
                 OR
           A reference to the classification array. In the latter case
           if there is a second argument and it evaluates to true,
           names will not be validated. NB: in any case, names are never
           validated anyway.

=cut

sub classification {
    my ($self, @vals) = @_;

    if (@vals) {
        if (ref($vals[0]) eq 'ARRAY') {
            @vals = @{$vals[0]};
        }
        
        # make sure the lineage contains us as first or second element
        # (lineage may have subspeces, species, genus ...)
        my $name = $self->node_name;
        if ($name && ($name ne $vals[0] && $name ne $vals[1]) &&
			       $name ne "$vals[1] $vals[0]") {
            $self->throw("The supplied lineage does not start near '$name'");
        }
        
        # create a lineage for ourselves
        my $db = Bio::DB::Taxonomy->new(-source => 'list', -names => [reverse @vals]);
        unless ($self->scientific_name) {
            # assume we're supposed to be the leaf of the supplied lineage
            $self->scientific_name($vals[0]);
        }
        unless ($self->rank) {
            # and that we are rank species
            $self->rank('species');
        }
        
        $self->db_handle($db);

        $self->{tree} = Bio::Tree::Tree->new(-node => $self);
        # some things want to freeze/thaw Bio::Species objects, but tree's
        # _root_cleanup_methods contains a CODE ref, delete it.
        #*** even if we don't delete the cleanup methods, we still get memory
        #    leak-like symtoms, and the actual cleanup causes a mass of
        #    warnings... needs investigation!
        delete $self->{tree}->{_root_cleanup_methods};
    }
    
    @vals = ();
    foreach my $node ($self->{tree}->get_lineage_nodes($self), $self) {
        unshift(@vals, $node->scientific_name || next);
    }
    weaken($self->{tree}->{'_rootnode'});
    return @vals;
}

=head2 ncbi_taxid

 Title   : ncbi_taxid
 Usage   : $obj->ncbi_taxid($newval)
 Function: Get/set the NCBI Taxon ID
 Returns : the NCBI Taxon ID as a string
 Args    : newvalue to set or undef to unset (optional)

=cut

=head2 common_name

 Title   : common_name
 Usage   : $self->common_name( $common_name );
           $common_name = $self->common_name();
 Function: Get or set the common name of the species
 Example : $self->common_name('human')
 Returns : The common name in a string
 Args    : String, which is the common name (optional)

=cut

=head2 division

 Title   : division
 Usage   : $obj->division($newval)
 Function: Genbank Division for a species
 Returns : value of division (a scalar)
 Args    : value of division (a scalar)

=cut

=head2 species

 Title   : species
 Usage   : $self->species( $species );
           $species = $self->species();
 Function: Get or set the scientific species name.
 Example : $self->species('Homo sapiens');
 Returns : Scientific species name as string
 Args    : Scientific species name as string

=cut

sub species {
    my ($self, $species) = @_;
    
	if ($species) {
		$self->{_species} = $species;
	}
	
	unless (defined $self->{_species}) {
		# work it out from our nodes
		my $species_taxon = $self->{tree}->find_node(-rank => 'species');
		unless ($species_taxon) {
			# just assume we are rank species
			$species_taxon = $self;
		}
		
		$species = $species_taxon->scientific_name;
		
		#
		# munge it like the Bio::SeqIO modules used to do
		# (more or less copy/pasted from old Bio::SeqIO::genbank, hence comments
		#  referring to 'ORGANISM' etc.)
		#
		
		my $root = $self->{tree}->get_root_node;
		unless ($root) {
            $self->{tree} = new Bio::Tree::Tree(-node => $species_taxon);
            delete $self->{tree}->{_root_cleanup_methods};
            $root = $self->{tree}->get_root_node;
        }
        
		my @spflds = split(' ', $species);
		if (@spflds > 1 && $root->node_name ne 'Viruses') {
			$species = undef;
			
			# does the next term start with uppercase?
			# yes: valid genus; no then unconventional
			# e.g. leaf litter basidiomycete sp. Collb2-39
			my $genus;
			if ($spflds[0] =~ m/^[A-Z]/) {
				$genus = shift(@spflds);
			}
			else {
				undef $genus;
			}
			
			my $sub_species;
			if (@spflds) {
				while (my $fld = shift @spflds) {
					$species .= "$fld ";
					# does it have subspecies or varieties?
					last if ($fld =~ m/(sp\.|var\.)/);
				}
				chop $species;	# last space
				$sub_species = join ' ',@spflds if(@spflds);
			}
			else {
				$species = 'sp.';
			}
			
			# does ORGANISM start with any words which make its genus undefined?
			# these are in @unkn_genus	
			# this in case species starts with uppercase so isn't caught above. 
			# alter common name if required
			my $unconv = 0; # is it unconventional species name?
			my @unkn_genus = ('unknown','unclassified','uncultured','unidentified');
			foreach (@unkn_genus) {
				if ($genus && $genus =~ m/$_/i)	{
					$species = $genus . " " . $species;
					undef $genus;
					$unconv = 1;
					last;
				}
				elsif ($species =~ m/$_/i)	{
					$unconv = 1;
					last;
				}
			}
			if (!$unconv && !$sub_species && $species =~ s/^(\w+)\s(\w+)$/$1/)	{
				# need to extract subspecies from conventional ORGANISM format.  
				# Will the 'word' in a two element species name
				# e.g. $species = 'thummi thummi' => $species='thummi' & 
				# $sub_species='thummi'
				$sub_species = $2;
			}
			
			$self->genus($genus) if $genus;
			$self->sub_species($sub_species) if $sub_species;
		}
		
		$self->{_species} = $species;
	}
	
	return $self->{_species};
}

=head2 genus

 Title   : genus
 Usage   : $self->genus( $genus );
           $genus = $self->genus();
 Function: Get or set the scientific genus name.
 Example : $self->genus('Homo');
 Returns : Scientific genus name as string
 Args    : Scientific genus name as string

=cut

sub genus {
    my ($self, $genus) = @_;
    
	if ($genus) {
        $self->{_genus} = $genus;
    }
	
	unless (defined $self->{_genus}) {
		my $genus_taxon = $self->{tree}->find_node(-rank => 'genus');
		unless ($genus_taxon) {
			# just assume our ancestor is rank genus
			$genus_taxon = $self->ancestor;
		}
		
		$self->{_genus} = $genus_taxon->scientific_name if $genus_taxon;
	}
	
	return $self->{_genus};
}

=head2 sub_species

 Title   : sub_species
 Usage   : $obj->sub_species($newval)
 Function: Get or set the scientific subspecies name.
 Returns : value of sub_species
 Args    : newvalue (optional)

=cut

sub sub_species {
    my ($self, $sub) = @_;
    
    unless (defined $self->{'_sub_species'}) {
        my $ss_taxon = $self->{tree}->find_node(-rank => 'subspecies');
        if ($ss_taxon) {
            if ($sub) {
                $ss_taxon->scientific_name($sub);
            }
            return $ss_taxon->scientific_name;
        }
    }
    
    # fall back to direct storage on self
    $self->{'_sub_species'} = $sub if $sub;
    return $self->{'_sub_species'};
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
    my ($self, $var) = @_;
    
    unless (defined $self->{'_variant'}) {
        my $var_taxon = $self->{tree}->find_node(-rank => 'variant');
        if ($var_taxon) {
            if ($var) {
                $var_taxon->scientific_name($var);
            }
            return $var_taxon->scientific_name;
        }
    }
    
    # fall back to direct storage on self
    $self->{'_variant'} = $var if $var;
    return $self->{'_variant'};
}

=head2 binomial

 Title   : binomial
 Usage   : $binomial = $self->binomial();
           $binomial = $self->binomial('FULL');
 Function: Returns a string "Genus species", or "Genus species subspecies",
           if the first argument is 'FULL' (and the species has a subspecies).
 Args    : Optionally the string 'FULL' to get the full name including
           the subspecies.

=cut

sub binomial {
    my ($self, $full) = @_;
    my $rank = $self->rank || 'no rank';
    
    my ($species, $genus) = ($self->species, $self->genus);
    unless (defined $species) {
        $species = 'sp.';
        $self->warn("requested binomial but classification was not set");
    }
    $genus = '' unless( defined $genus);
    
    $species =~ s/$genus\s+//;
    
    my $bi = "$genus $species";
    if (defined($full) && $full =~ /full/i) { 
        my $ssp = $self->sub_species;
        if ($ssp) {
            $ssp =~ s/$bi\s+//;
            $ssp =~ s/$species\s+//;
            $bi .= " $ssp";
        }
    }
    return $bi;
}

=head2 validate_species_name

 Title   : validate_species_name
 Usage   : $result = $self->validate_species_name($string);
 Function: Validate the species portion of the binomial
 Args    : string
 Notes   : The string following the "genus name" in the NCBI binomial
           is so variable that it's not clear that this is a useful
           function. Consider the binomials 
           "Simian 11 rotavirus (serotype 3 / strain SA11-Patton)",
           or "St. Thomas 3 rotavirus", straight from GenBank.
           This is particularly problematic in microbes and viruses.
           As such, this isn't actually used automatically by any Bio::Species
           method.
=cut

sub validate_species_name {
    my( $self, $string ) = @_;

    return 1 if $string eq "sp.";
	return 1 if $string =~ /strain/;
    return 1 if $string =~ /^[a-z][\w\s-]+$/i;
    $self->throw("Invalid species name '$string'");
}

sub validate_name {
    return 1;
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
    my($self) = shift;
    return $self->{'_organelle'} = shift if @_;
    return $self->{'_organelle'};
}

sub dont_DESTROY {
    my $self = shift;
    $self->{tree}->cleanup_tree if $self->{tree};
    delete $self->{tree};
    $self->node_cleanup;
}

1;
