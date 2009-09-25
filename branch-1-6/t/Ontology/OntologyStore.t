# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 6,
			   -requires_module => 'Graph',
			   -requires_networking => 1);
	
	use_ok('Bio::Ontology::OntologyStore');
}

ok my $store = Bio::Ontology::OntologyStore->get_instance;
SKIP: {
	my $ontology;
	eval {$ontology = $store->get_ontology(-name => 'Sequence Ontology');};
	skip "Couldn't get sequence ontology, network problems? Skipping these tests", 4 if $@;
	ok(scalar($ontology->get_root_terms()) == 1);
	my($txt) = $ontology->find_terms(-name => 'transcript');
	is $txt->identifier, 'SO:0000673';
	is $txt->name, 'transcript';
	is $txt->definition, 'An RNA synthesized on a DNA or RNA template by an RNA polymerase.';
}
