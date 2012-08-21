use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    test_begin( -tests => 42 );
    use_ok('Bio::DB::Taxonomy');
    use_ok('Bio::Tree::Tree');
}


my ($db, $id, @ids, $node, $node2, $ancestor, @children, %common_names, $tree, @descs);


# Test Bio::DB::Taxonomy::silva

ok $db = Bio::DB::Taxonomy->new( -source => 'silva' );

isa_ok $db, 'Bio::DB::Taxonomy::silva';
isa_ok $db, 'Bio::DB::Taxonomy::list';
isa_ok $db, 'Bio::DB::Taxonomy';

ok $db = Bio::DB::Taxonomy->new(
   -source   => 'silva',
   -taxofile => test_input_file('taxonomy', 'silva_SSURef_108_tax_silva_trunc.fasta'),
);

@ids = $db->get_taxonid('Homo sapiens');
is scalar @ids, 0;

@ids = $db->get_taxonid('Rattus norvegicus');
is scalar @ids, 1;
$id = $ids[0];

ok $node = $db->get_taxon($id);
is $node->id, $id;
is $node->object_id, $node->id;
is $node->ncbi_taxid, undef;
is $node->rank, undef;
is $node->parent_id, 'sv72';
is $node->node_name, 'Rattus norvegicus';
is $node->scientific_name, $node->node_name;

is ${$node->name('scientific')}[0], $node->node_name;

# it has a common name in Silva, but Bio::DB::Taxonomy::silva does not record it
%common_names = map { $_ => 1 } $node->common_names;
is scalar keys %common_names, 0;

is $node->division, undef;
is $node->genetic_code, undef;
is $node->mitochondrial_genetic_code, undef;

# briefly test some Bio::Tree::NodeI methods
ok $ancestor = $node->ancestor;
is $ancestor->scientific_name, '';
ok $ancestor = $ancestor->ancestor;
is $ancestor->scientific_name, 'Rattus';
ok $ancestor = $ancestor->ancestor;
is $ancestor->scientific_name, 'Murinae';

# unless set explicitly, Bio::Taxon doesn't return anything for
# each_Descendent; must ask the database directly
@ids = $db->get_taxonid('Metazoa');
is scalar @ids, 1;
$id = $ids[0];
ok $node = $db->get_taxon($id);
ok @children = $node->db_handle->each_Descendent($node);
is scalar @children, 3; # Chordata, Platyhelminthes, Metazoa

# do some trickier things...
ok $node2 = $db->get_taxon('sv112');
is $node2->scientific_name, 'Mucoromycotina';

# briefly check that we can use some Tree methods
$tree = Bio::Tree::Tree->new();
is $tree->get_lca($node, $node2)->scientific_name, 'Eukaryota';

# can we actually form a Tree and use other Tree methods?
ok $tree = Bio::Tree::Tree->new(-node => $node2);
is $tree->number_nodes, 5;
is $tree->get_nodes, 5;

# check that getting the ancestor still works now we have explitly set the
# ancestor by making a Tree
is $node2->ancestor->scientific_name, 'Basal fungal lineages';

# we can recursively fetch all descendents of a taxon
my $lca = $db->get_taxon( -name => 'Liliopsida' );
ok @descs = $db->each_Descendent($lca);
is scalar @descs, 1;
@descs = $db->get_all_Descendents($lca);
is scalar @descs, 9;

