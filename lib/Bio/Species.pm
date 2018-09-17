#
# BioPerl module for Bio::Species
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by James Gilbert <jgrg@sanger.ac.uk>
# Reimplemented by Sendu Bala <bix@sendu.me.uk>
# Re-reimplemented by Chris Fields <cjfields - at - bioperl dot org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Species - Generic species object.  

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

B<NOTE: This class is planned for deprecation in favor of the simpler Bio::Taxon.
Please use that class instead.>

Provides a very simple object for storing phylogenetic information. The
classification is stored in an array, which is a list of nodes in a phylogenetic
tree. Access to getting and setting species and genus is provided, but not to
any of the other node types (eg: "phylum", "class", "order", "family"). There's
plenty of scope for making the model more sophisticated, if this is ever needed.

A methods are also provided for storing common names, and subspecies.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

James Gilbert email B<jgrg@sanger.ac.uk>

=head1 CONTRIBUTORS

Sendu Bala, bix@sendu.me.uk
Chris Fields, cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' Let the code begin...

package Bio::Species;
use strict;
use warnings;

use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::Taxon;
use base qw(Bio::Root::Root Bio::Tree::NodeI);

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
    
    # Bio::Species is now just a proxy object that just observes the NodeI
    # interface methods but delegates them to the proper classes (Bio::Taxon and
    # Bio::Tree::Tree). This will be surplanted by the much simpler
    # Bio::Taxon/Bio::DB::Taxonomy modules in the future.
    
    # Using a proxy allows proper GC w/o using weaken().  This just wraps the
    # older instances, which have no reciprocal refs (thus no circular refs).
    # This can then run proper cleanup
    
    $self->taxon(Bio::Taxon->new(@args));
    
    my ($org, $sp, $var, $classification) =
        $self->_rearrange([qw(ORGANELLE
                            SUB_SPECIES
                            VARIANT
                            CLASSIFICATION)], @args);
    
    if (defined $classification && ref($classification) eq "ARRAY" && @{$classification}) {
        $self->classification(@$classification);
    }
    else {
        $self->tree(Bio::Tree::Tree->new());
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

    my $taxon = $self->taxon;

    if (@vals) {
        if (ref($vals[0]) eq 'ARRAY') {
            @vals = @{$vals[0]};
        }
        
        $vals[1] ||= '';
        # make sure the lineage contains us as first or second element
        # (lineage may have subspecies, species, genus ...)
        my $name = $taxon->node_name;
        my ($genus, $species) = (quotemeta($vals[1]), quotemeta($vals[0]));
        if ($name && 
           ($name !~ m{$species}i && $name !~ m{$genus}i) && 
            $name !~ m{$genus $species}i) {
            if ($name =~ /^$genus $species\s*(.+)/) {
                # just assume the problem is someone tried to make a Bio::Species starting at subspecies
                #*** no idea if this is appropriate! just a possible fix related to bug 2092
                $self->sub_species($1);
                $name = $taxon->node_name("$vals[1] $vals[0]");
            }
            else {
                $self->warn("The supplied lineage does not start near '$name' (I was supplied '".join(" | ", @vals)."')");
            }
        }
        
        # create a lineage for ourselves
        my $db = Bio::DB::Taxonomy->new(-source => 'list', -names => [reverse @vals]);
        unless ($taxon->scientific_name) {
            # assume we're supposed to be the leaf of the supplied lineage
            $self->taxon->scientific_name($vals[0]);
        }
        unless ($taxon->rank) {
            # and that we are rank species
            $taxon->rank('species');
        }
        
        $taxon->db_handle($db);
        
        $self->tree(Bio::Tree::Tree->new(-node => $taxon));
    }
    
    @vals = ();
    foreach my $node ($self->tree->get_lineage_nodes($taxon), $taxon) {
        unshift(@vals, $node->scientific_name || next);
    }
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
 Function: Get or set the species name.
           Note that this is  NOT genus and species
           -- use $self->binomial() for that.
 Example : $self->species('sapiens');
 Returns : species name as string (NOT genus and species)
 Args    : species name as string (NOT genus and species)

=cut

sub species {
    my ($self, $species) = @_;
    
	if ($species) {
		$self->{_species} = $species;
	}

	unless (defined $self->{_species}) {
		# work it out from our nodes
		my $species_taxon = $self->tree->find_node(-rank => 'species');
		unless ($species_taxon) {
			# just assume we are rank species
			$species_taxon = $self->taxon;
		}

		$species = $species_taxon->scientific_name;
        
		#
		# munge it like the Bio::SeqIO modules used to do
		# (more or less copy/pasted from old Bio::SeqIO::genbank, hence comments
		#  referring to 'ORGANISM' etc.)
		#

		my $root = $self->tree->get_root_node;
		unless ($root) {
            $self->tree(Bio::Tree::Tree->new(-node => $species_taxon));
            $root = $self->tree->get_root_node;
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

    # TODO: instead of caching the raw name, cache the actual node instance.
    if ($genus) {
        $self->{_genus} = $genus;
    }
    unless (defined $self->{_genus}) {
        my $genus_taxon = $self->tree->find_node(-rank => 'genus');
        unless ($genus_taxon) {
            # just assume our ancestor is rank genus
            $genus_taxon = $self->taxon->ancestor;
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
    
    # TODO: instead of caching the raw name, cache the actual node instance.
    if (!defined $self->{'_sub_species'}) {
        my $ss_taxon = $self->tree->find_node(-rank => 'subspecies');
        if ($ss_taxon) {
            if ($sub) {
                $ss_taxon->scientific_name($sub);
                
                # *** weakening ref to our root node in species() to solve a
                # memory leak means that we have a subspecies taxon to set
                # during the first call to species(), but it has vanished by
                # the time a user subsequently calls sub_species() to get the
                # value. So we 'cheat' and just store the subspecies name in
                # our self hash, instead of the tree. Is this a problem for
                # a Species object? Can't decide --sendu
                
                # This can now be changed to deal with this information on the
                # fly.  For now, the caching remains, but maybe we should just
                # let these things deal with mutable data as needed? -- cjfields
                
                $self->{'_sub_species'} = $sub;
            }
            return $ss_taxon->scientific_name;
        }
        else {
            # should we create a node here to be added to the tree?
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
    
    # TODO: instead of caching the raw name, cache the actual node instance.
    if (!defined $self->{'_variant'}) {
        my $var_taxon = $self->tree->find_node(-rank => 'variant');
        if ($var_taxon) {
            if ($var) {
                $var_taxon->scientific_name($var);
            }
            return $var_taxon->scientific_name;
        }
        else {
            # should we create a node here to be added to the tree?
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
 Note    : This is just munged from the taxon() name

=cut

sub binomial {
    my ($self, $full) = @_;
    my $rank = $self->taxon->rank || 'no rank';
    
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
 Notes   : The string following the "genus name" in the NCBI binomial is so
           variable that it's not clear that this is a useful function. Consider
           the binomials "Simian 11 rotavirus (serotype 3 / strain
           SA11-Patton)", or "St. Thomas 3 rotavirus", straight from GenBank.
           This is particularly problematic in microbes and viruses. As such,
           this isn't actually used automatically by any Bio::Species method.

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
 Note    : TODO: We currently do not know where the organelle definition will
           eventually go.  This is stored in the source seqfeature, though,
           so the information isn't lost.

=cut

sub organelle {
    my($self) = shift;
    return $self->{'_organelle'} = shift if @_;
    return $self->{'_organelle'};
}

=head2 Delegation

The following methods delegate to the internal Bio::Taxon instance. This is
mainly to allow code continue using older methods, with the mind to migrate to
using Bio::Taxon and related methods when this class is deprecated.

=cut 

sub node_name {shift->taxon->node_name(@_)}
sub scientific_name {shift->taxon->node_name(@_)}

sub id {shift->taxon->id(@_)}
sub object_id {shift->taxon->id(@_)}
sub ncbi_taxid {shift->taxon->ncbi_taxid(@_)}
sub rank {shift->taxon->rank(@_)}
sub division {shift->taxon->division(@_)}

sub common_names {shift->taxon->common_names(@_)}
sub common_name {shift->taxon->common_names(@_)}

sub genetic_code {shift->taxon->genetic_code(@_)}
sub mitochondrial_genetic_code {shift->taxon->mitochondrial_genetic_code(@_)}

sub create_date { shift->taxon->create_date(@_)}
sub pub_date { shift->taxon->pub_date(@_)}
sub update_date { shift->taxon->update_date(@_)}

sub db_handle { shift->taxon->db_handle(@_)}

sub parent_id { shift->taxon->parent_id(@_)}
sub parent_taxon_id { shift->taxon->parent_id(@_)}

sub version { shift->taxon->version(@_)}
sub authority { shift->taxon->authority(@_)}
sub namespace { shift->taxon->namespace(@_)}

sub ancestor { shift->taxon->ancestor(@_)}
sub get_Parent_Node { shift->taxon->get_Parent_Node(@_)}
sub each_Descendent { shift->taxon->each_Descendent(@_)}
sub get_Children_Nodes { shift->taxon->get_Children_Nodes(@_)}
sub remove_Descendant { shift->taxon->remove_Descendant(@_)}

sub name { shift->taxon->name(@_)}

=head2 taxon

 Title    : taxon
 Usage    : $obj->taxon
 Function : retrieve the internal Bio::Taxon instance
 Returns  : A Bio::Taxon. If one is not previously set,
            an instance is created lazily
 Args     : Bio::Taxon (optional)
 
=cut

sub taxon {
    my ($self, $taxon) = @_;
    if (!$self->{taxon} || $taxon) {
        $taxon ||= Bio::Taxon->new();
        $self->{taxon} = $taxon;
    }
    $self->{taxon};
}

=head2 tree

 Title    : tree
 Usage    : $obj->tree
 Function : Returns a Bio::Tree::Tree object
 Returns  : A Bio::Tree::Tree. If one is not previously set,
            an instance is created lazily
 Args     : Bio::Tree::Tree (optional)
 
=cut

sub tree {
    my ($self, $tree) = @_;
    if (!$self->{tree} || $tree) {
        $tree ||= Bio::Tree::Tree->new();
        delete $tree->{_root_cleanup_methods};
        $self->{tree} = $tree;
    }
    $self->{tree};
}

sub DESTROY {
    my $self = shift;
    $self->tree->cleanup_tree;
    delete $self->{tree};
    $self->taxon->node_cleanup;
}

1;
