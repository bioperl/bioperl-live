# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 55,
			   -requires_modules => [qw(XML::Parser::PerlSAX
									    XML::Parser
										Graph::Directed)]);
	
	use_ok('Bio::OntologyIO');
}

my $ipp = Bio::OntologyIO->new( -format => 'interpro',
										  -file => test_input_file('interpro_short.xml'),
										  -ontology_engine => 'simple' );
isa_ok ($ipp, 'Bio::OntologyIO::InterProParser');

my $ip;
while(my $ont = $ipp->next_ontology()) {
    # there should be only one ontology
    is ($ip, undef);
    $ip = $ont;
}
#print $ip->to_string."\n";
my @rt = sort { $a->name cmp $b->name; } $ip->get_root_terms();

# there should be 8 root terms: InterPro Domain, InterPro Family,
# InterPro Repeat, and InterPro PTM (Post Translational Modification),
# Active_site, Binding_site, Conserved_site and Region

is (scalar(@rt), 8, 'There are 8 root InterPro terms');

# every InterPro term should have an ontology,
foreach ($ip->get_leaf_terms, @rt) {
	isa_ok ($_->ontology, 'Bio::Ontology::Ontology');
	is ($_->ontology->name, "InterPro",
		 "term ".$_->name." in ontology InterPro");
}

# there are 10 fully instantiated InterPro terms in total, which should be returned as the leafs
is (scalar($ip->get_leaf_terms()), 10);
# roots and leafs together:
is (scalar($ip->get_all_terms()), 15);

# descendants and children (synonymous here because of depth 1)
# note that the sort should have placed Domain first and Family second
is (scalar($ip->get_descendant_terms($rt[3])), 4); # 4 InterPro Domains
is (scalar($ip->get_child_terms($rt[3])), 4);      # dto.
is (scalar($ip->get_descendant_terms($rt[4])), 3); # 3 Interpro Family
is (scalar($ip->get_child_terms($rt[4])), 3);      # dto.

# test for ancestors and parents (synonymous here because of depth 1)
foreach my $t ($ip->get_leaf_terms) {
	# every InterPro term has exactly one parent - namely either 
	# Domain, Region, Family, Repeat, or PTM(Post Transl. Modification)

	if (!($t->identifier eq "Repeat" || $t->identifier eq "PTM" || $t->identifier eq "Region"
			|| $t->identifier =~ '_site' )) {
		is (scalar($ip->get_parent_terms($t)), 1, $t->name . " term has one parent");
		is (scalar($ip->get_ancestor_terms($t)), 1, $t->name . " term has one ancestor");
	}
}

# test for secondary accession map
is(scalar(keys %{$ipp->secondary_accessions_map}), 2, 'secondary accession map has 2 keys');
