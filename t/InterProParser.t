# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 47,
			   -requires_modules => [qw(XML::Parser::PerlSAX
									    XML::Parser
										Graph::Directed)]);
	
	use_ok('Bio::OntologyIO');
}

my $ipp = Bio::OntologyIO->new( -format => 'interpro',
										  -file => test_input_file('interpro_short.xml'),
										  -ontology_engine => 'simple' );
ok ($ipp);

my $ip;
while(my $ont = $ipp->next_ontology()) {
    # there should be only one ontology
    is ($ip, undef);
    $ip = $ont;
}
#print $ip->to_string."\n";
my @rt = sort { $a->name cmp $b->name; } $ip->get_root_terms();

# there should be 4 root terms: InterPro Domain, InterPro Family,
# InterPro Repeat, and InterPro PTM (Post Translational Modification)
#
# I added 2 more terms, Active_site and Binding_site. -Juguang
is (scalar(@rt), 6);

# every InterPro term should have an ontology,
foreach ($ip->get_leaf_terms, @rt) {
	ok ($_->ontology);
	is ($_->ontology->name, "InterPro",
		 "term ".$_->name." in ontology InterPro");
}

# there are 6 fully instantiated InterPro terms in total, which should be returned as the leafs
is (scalar($ip->get_leaf_terms()), 8);
# roots and leafs together:
is (scalar($ip->get_all_terms()), 13);

# descendants and children (synonymous here because of depth 1)
# note that the sort should have placed Domain first and Family second
is (scalar($ip->get_descendant_terms($rt[2])), 4); # 4 InterPro Domains
is (scalar($ip->get_child_terms($rt[2])), 4);      # dto.
is (scalar($ip->get_descendant_terms($rt[3])), 3); # 3 Interpro Family
is (scalar($ip->get_child_terms($rt[3])), 3);      # dto.

# test for ancestors and parents (synonymous here because of depth 1)
foreach my $t ($ip->get_leaf_terms) {
	# every InterPro term has exactly one parent - namely either 
	# Domain, Family, Repeat, or PTM(Post Transl. Modification)
	if (!($t->identifier eq "Repeat" || $t->identifier eq "PTM" 
			|| $t->identifier eq 'Active_site' || $t->identifier eq 'Binding_site')) {
		is (scalar($ip->get_parent_terms($t)), 1);
		is (scalar($ip->get_ancestor_terms($t)), 1);
	}
}

# test for secondary accession map
is(scalar(keys %{$ipp->secondary_accessions_map}), 2);
