# -*-Perl-*- Test Harness script for Bioperl
# $Id$

# These modules are now deprecated, don't bother testing them. --sendubala

## I am pretty sure this module is going the way of the dodo bird so 
## I am not sure how much work to put into fixing the tests/module
## --jasonstajich

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 0);
}

if (0) {
	use Bio::Taxonomy::Taxon;
	ok(1);
	
	
	ok my $taxonL = Bio::Taxonomy::Taxon->new;
	ok $taxonL->description('this could be anything');
	ok $taxonL->taxon('could this be called name?');
	ok $taxonL->id('could this be called taxid?');
	skip 1, $taxonL->branch_length('should accept only numerical values?');
	ok  $taxonL->branch_length(5);
	
	ok $taxonL->id('could this be called taxid?');
	ok $taxonL->rank('species');
	ok $taxonL->rank, 'species';
	# ok $taxonL->has_rank, 'species'; #why two methods that do mostly the same thing, but work differently?
	
	skip 1, $taxonL->rank('foo is not a rank, class variable @RANK not initialised'); 
	ok $taxonL->to_string, '"could this be called taxid?":5';
	
	my $taxonR = Bio::Taxonomy::Taxon->new();
	
	my $taxon = Bio::Taxonomy::Taxon->new(-id =>'ancient', -taxon => 'genus');
	ok $taxon->id(), 'ancient'; 
	ok $taxon->taxon(), 'genus'; 
	ok $taxon->internal_id, 2;
	ok $taxonL->internal_id, 0; # would not it be better to start numebering from 1?
	ok $taxon->add_Descendent($taxonL);
	$taxon->add_Descendent($taxonR);
	
	ok  scalar $taxon->each_Descendent, 2;  # dies
	ok $taxon->remove_Descendent($taxonR); # better to return number of Descendants removed
	
	ok $taxon->remove_all_Descendents();
	
	
	$taxon->add_Descendent($taxonL);
	ok $taxonL->ancestor->id, 'ancient';
	ok $taxonL->branch_length(5);
	
	
	ok $taxonL->is_Leaf, 1;
	ok $taxon->is_Leaf, 0;
	ok $taxon->height, 6;
	ok $taxonL->height, 5;
	ok $taxon->invalidate_height, undef;
	ok $taxonL->classify(1), 2;
	skip(1,"skip classify weirdness");
	# ok $taxonL->classify(0), 2, 'ancestor has rank, but implementation prevents showing anything more than one value';
	skip(1,"skip classify weirdness");
	#ok $taxonL->has_rank, 1, 'documentation claims this returns a boolean; and that it queries ancestors rank?, needs an agrument but does not test it';
	skip(1,"skip classify weirdness");
	#ok $taxonL->has_rank('species'), 1;
	
	#ok $taxon->has_taxon(); # why docs and code talk about ancestor?
	#ok $taxonL->has_taxon('genus');  returns undef or oan object, not boolean
	
	ok $taxon->distance_to_root, 0;
	ok $taxonL->distance_to_root, 1;
	#ok $taxonL->recent_common_ancestor($taxon)->id, 'ancient';
	
	
	
	#use Data::Dumper;
	#print Dumper  $taxonL->classify();
	skip(1, 'Skip this weird function');
	# ok $taxonL->has_rank('species'), 1;
	#ok my $species = $taxonL->species;
	
	
	
	
	
	##################################################################################################
	
	# tests for Bio::Taxonomy::Tree;
	# code from synopsis
	
	use Bio::Species;
	use Bio::Taxonomy::Tree;
	use Bio::Taxonomy;
	
	my $human=Bio::Species->new();
	my $chimp=Bio::Species->new();
	my $bonobo=Bio::Species->new();
	
	$human->classification(qw( sapiens Homo Hominidae
							   Catarrhini Primates Eutheria
							   Mammalia Euteleostomi Vertebrata 
							   Craniata Chordata
							   Metazoa Eukaryota ));
	$chimp->classification(qw( troglodytes Pan Hominidae
							   Catarrhini Primates Eutheria
							   Mammalia Euteleostomi Vertebrata 
							   Craniata Chordata
							   Metazoa Eukaryota ));
	$bonobo->classification(qw( paniscus Pan Hominidae
								Catarrhini Primates Eutheria
								Mammalia Euteleostomi Vertebrata 
								Craniata Chordata
								Metazoa Eukaryota ));
	
	# ranks passed to $taxonomy match ranks of species
	my @ranks = ('superkingdom','kingdom','phylum','subphylum',
				 'no rank 1','no rank 2','class','no rank 3','order',
				 'suborder','family','genus','species');
	
	my $taxonomy=Bio::Taxonomy->new(-ranks => \@ranks,
								   -method => 'trust',
								   -order => -1);
	
	
	ok my $tree1=Bio::Taxonomy::Tree->new();
	my $tree2=Bio::Taxonomy::Tree->new();
	
	$tree1->make_species_branch($human,$taxonomy);
	$tree2->make_species_branch($chimp,$taxonomy);
	
	my ($homo_sapiens) = $tree1->get_leaves;
	ok ref $homo_sapiens, 'Bio::Taxonomy::Taxon';
	
	ok $tree1->splice($tree2);
	
	ok $tree1->add_species($bonobo,$taxonomy);
	
	
	ok join (", ", map {$_->taxon} $tree1->get_leaves), 'Homo sapiens, Pan troglodytes, Pan paniscus';
	ok $tree1->remove_branch($homo_sapiens);
	ok join (", ", map {$_->taxon} $tree1->get_leaves), 'Pan troglodytes, Pan paniscus';
}
