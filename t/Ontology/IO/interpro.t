# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 69,
               -requires_modules => [qw(XML::Parser::PerlSAX
                                        XML::Parser
                                        Graph::Directed)]);
    use_ok('Bio::OntologyIO');
}
use Data::Dumper;

my $ipp = Bio::OntologyIO->new( -format          => 'interpro',
                                -file            => test_input_file('interpro_short.xml'),
                                -ontology_engine => 'simple' );
isa_ok ($ipp, 'Bio::OntologyIO::InterProParser');

my $ip;
while(my $ont = $ipp->next_ontology) {
    # there should be only one ontology
    is ($ip, undef);
    $ip = $ont;
}
# we grep for defined values because we don't want a list of undefined values to be pass any test
my @leaves = $ip->get_leaf_terms;
ok( ( grep { defined } map { $_->get_dbxrefs } @leaves), 'get_dbxrefs on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('member_list') }       @leaves), 'get_dbxrefs(member_list) on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('sec_list') }          @leaves), 'get_dbxrefs(sec_list) on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('class_list') }        @leaves), 'get_dbxrefs(class_list) on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('pub_list') }          @leaves), 'get_dbxrefs(pub_list) on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('example_list') }      @leaves), 'get_dbxrefs(example_list) on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_dbxrefs('external_doc_list') } @leaves), 'get_dbxrefs(external_doc_list) on leaf terms is non-empty');

ok( ( grep { defined } map { $_->get_members }            @leaves), 'get_members on leaf terms is non-empty');
ok( ( grep { defined } map { $_->class_list }             @leaves), 'class_list on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_examples }           @leaves), 'get_examples on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_external_documents } @leaves), 'get_external_documents on leaf terms is non-empty');
ok( ( grep { defined } map { $_->get_references }         @leaves), 'get_references on leaf terms is non-empty');
ok( ( grep { defined } map { $_->protein_count }          @leaves), 'protein_count on leaf terms is non-empty');

# this could greatly be improved (note: @leaves elements currently come in random order)
like( $leaves[0]->to_string, qr/-- InterPro id:/, 'to_string looks reasonable');

# there should be 8 root terms: InterPro Domain, InterPro Family,
# InterPro Repeat, and InterPro PTM (Post Translational Modification),
# Active_site, Binding_site, Conserved_site and Region

my @rt = sort { $a->name cmp $b->name; } $ip->get_root_terms();
is (scalar(@rt), 8, 'There are 8 root InterPro terms');

# every InterPro term should have an ontology,
foreach my $term ($ip->get_leaf_terms, @rt) {
    isa_ok ($term->ontology, 'Bio::Ontology::Ontology');
    is ($term->ontology->name, "InterPro", "term ".$term->name." in ontology InterPro");
}

# there are 10 fully instantiated InterPro terms in total, which should be returned as the leafs
is (scalar($ip->get_leaf_terms), 10);
# roots and leafs together:
is (scalar($ip->get_all_terms), 15);

# descendants and children (synonymous here because of depth 1)
# note that the sort should have placed Domain first and Family second
is (scalar($ip->get_descendant_terms($rt[3])), 4); # 4 InterPro Domains
is (scalar($ip->get_child_terms($rt[3])),      4); # dto.
is (scalar($ip->get_descendant_terms($rt[4])), 3); # 3 Interpro Family
is (scalar($ip->get_child_terms($rt[4])),      3); # dto.

# test for ancestors and parents (synonymous here because of depth 1)
foreach my $t ($ip->get_leaf_terms) {
    # every InterPro term has exactly one parent - namely either
    # Domain, Region, Family, Repeat, or PTM(Post Transl. Modification)

    if (not (   $t->identifier eq "Repeat" or $t->identifier eq "PTM"
             or $t->identifier eq "Region" or $t->identifier =~ '_site' )
        ) {
        is (scalar($ip->get_parent_terms($t)),   1, $t->name . " term has one parent");
        is (scalar($ip->get_ancestor_terms($t)), 1, $t->name . " term has one ancestor");
    }
}

# test for secondary accession map
is(scalar(keys %{$ipp->secondary_accessions_map}), 2, 'secondary accession map has 2 keys');
