use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin( -tests => 103 );

    use_ok('Bio::DB::Taxonomy');
    use_ok('Bio::Tree::Tree');
}


my ($db, $id, @ids, $node, $node2, $ancestor, @children, %common_names, $tree, @descs);


# Test Bio::DB::Taxonomy::greengenes

ok $db = Bio::DB::Taxonomy->new( -source => 'greengenes' );

isa_ok $db, 'Bio::DB::Taxonomy::greengenes';
isa_ok $db, 'Bio::DB::Taxonomy::list';
isa_ok $db, 'Bio::DB::Taxonomy';

ok $db = Bio::DB::Taxonomy->new(
   -source   => 'greengenes',
   -taxofile => test_input_file('taxonomy', 'greengenes_taxonomy_16S_candiv_gg_2011_1.txt'),
);


####
use Data::Dumper;
print Dumper($db);
####

@ids = $db->get_taxonid('s__Bacteroides uniformis');
is scalar @ids, 1;
$id = $ids[0];

@ids = $db->get_taxonid('Homo sapiens');
is scalar @ids, 0;

ok $node = $db->get_taxon($id);
is $node->id, $id;
is $node->object_id, $node->id;

is $node->ncbi_taxid, undef;
### record otu_id??

is $node->rank, 'species';

###print "parent id: ".$node->parent_id."\n";
###is $node->parent_id, 9605;

        
is $node->node_name, 's__Bacteroides uniformis';

is $node->scientific_name, $node->node_name;
### remove 'x__' from scientific names?

is ${$node->name('scientific')}[0], $node->node_name;

%common_names = map { $_ => 1 } $node->common_names;
is scalar keys %common_names, 0;

is $node->division, undef;
is $node->genetic_code, undef;
is $node->mitochondrial_genetic_code, undef;

# briefly test some Bio::Tree::NodeI methods
ok $ancestor = $node->ancestor;
is $ancestor->scientific_name, 'g__Bacteroides';

# unless set explicitly, Bio::Taxon doesn't return anything for
# each_Descendent; must ask the database directly
ok @children = $ancestor->db_handle->each_Descendent($ancestor);
is scalar @children, 2;

# do some trickier things...
ok $node2 = $db->get_taxon('list112');
is $node2->scientific_name, 'o__Synergistales';

# briefly check that we can use some Tree methods
$tree = Bio::Tree::Tree->new();
is $tree->get_lca($node, $node2)->scientific_name, 'k__Bacteria';

# can we actually form a Tree and use other Tree methods?
ok $tree = Bio::Tree::Tree->new(-node => $node);
#is $tree->number_nodes, 30;
#is $tree->get_nodes, 30;
is $tree->find_node(-rank => 'genus')->scientific_name, 'g__Bacteroides';

# check that getting the ancestor still works now we have explitly set the
# ancestor by making a Tree
is $node->ancestor->scientific_name, 'g__Bacteroides';


## entrez isn't as good at searching as flatfile, so we have to special-case
#my @ids = $db->get_taxonids('Chloroflexi');
#$db eq $db_entrez ? (is @ids, 1) : (is @ids, 2);
#$id = $db->get_taxonids('Chloroflexi (class)');
#is $id, 32061;
#
#@ids = $db->get_taxonids('Rhodotorula');
#cmp_ok @ids, '>=' , 8;
#@ids = $db->get_taxonids('Rhodotorula <Microbotryomycetidae>');
#is @ids, 1;
#is $ids[0], 231509;


# we can recursively fetch all descendents of a taxon

my $lca = $db->get_taxon( -name => 'f__Enterobacteriaceae' );
ok @descs = $db->each_Descendent($lca);
is scalar @descs, 2;
@descs = $db->get_all_Descendents($lca);
is scalar @descs, 3;


