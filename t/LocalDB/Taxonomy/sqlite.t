# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    #test_begin(
    #    
    #    -requires_modules => [qw(DBI DBD::SQLite )]
    #);

    use_ok('Bio::DB::Taxonomy');
    use_ok('Bio::Tree::Tree');
}

my $temp_dir = test_output_dir();

# TODO: run basic tests making sure that a database is not regenerated if
# present or unless forced

ok my $db_flatfile = Bio::DB::Taxonomy->new(
    -source    => 'sqlite',
    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
    -namesfile => test_input_file('taxdump', 'names.dmp'),
);
isa_ok $db_flatfile, 'Bio::DB::Taxonomy::sqlite';
isa_ok $db_flatfile, 'Bio::DB::Taxonomy';

# By not specifying a '-directory' argument, index files go to a temporary
# folder ($Bio::Root::IO::TEMPDIR, such as 'C:\Users\USER\AppData\Local\Temp'),
# and are implied to be temporary. So test the ability of flatfile->DESTROY to
# remove the temporary index files at object destruction (this also affects files
# in "test_output_dir()", since the folder is created inside the temporary folder)
#no warnings qw(once); # silence 'Name "$Bio::Root::IO::TEMPDIR" used only once'
#
#is $db_flatfile->{index_directory}, $Bio::Root::IO::TEMPDIR, 'removal of temporary index files: no -directory';

#$db_flatfile->DESTROY;
#ok not -e ($db_flatfile->{index_directory} . '/id2names');
#ok not -e ($db_flatfile->{index_directory} . '/names2id');
#ok not -e ($db_flatfile->{index_directory} . '/nodes');
#ok not -e ($db_flatfile->{index_directory} . '/parents');

## Test removal of temporary index files from test_output_dir folder
## (since test_output_dir() =~ m/^$Bio::Root::IO::TEMPDIR/)
#ok $db_flatfile = Bio::DB::Taxonomy->new(
#    -source    => 'flatfile',
#    -directory => $temp_dir,
#    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
#    -namesfile => test_input_file('taxdump', 'names.dmp'),
#    -force     => 1,
#);
#is $db_flatfile->{index_directory}, $temp_dir, 'removal of temporary index files: test_output_dir()';
#$db_flatfile->DESTROY;
#ok not -e ($db_flatfile->{index_directory} . '/id2names');
#ok not -e ($db_flatfile->{index_directory} . '/names2id');
#ok not -e ($db_flatfile->{index_directory} . '/nodes');
#ok not -e ($db_flatfile->{index_directory} . '/parents');
#
# Generate the object (and the files) again for the remaining tests

ok my $db = Bio::DB::Taxonomy->new(
    -source    => 'sqlite',
    -directory => $temp_dir,
    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
    -namesfile => test_input_file('taxdump', 'names.dmp'),
    -force     => 1,
);

my $id;

# taxid data in the nodes.dmp file should be unique, we ignore repeated values
# if seen

is $db->get_num_taxa, 188;

lives_ok {$id = $db->get_taxonid('Homo sapiens')};

is $id, 9606;

## easy test on human, try out the main Taxon methods
#ok $n = $db->get_taxon(9606);
#is $n->id, 9606;
#is $n->object_id, $n->id;
#is $n->ncbi_taxid, $n->id;
#is $n->parent_id, 9605;
#is $n->rank, 'species';
#
#is $n->node_name, 'Homo sapiens';
#is $n->scientific_name, $n->node_name;
#is ${$n->name('scientific')}[0], $n->node_name;
#
#my %common_names = map { $_ => 1 } $n->common_names;
#is keys %common_names, 3, ref($db).": common names";
#ok exists $common_names{human};
#ok exists $common_names{man};
#
#is $n->division, 'Primates';
#is $n->genetic_code, 1;
#is $n->mitochondrial_genetic_code, 2;
## these are entrez-only, data not available in dmp files
#if ($db eq $db_entrez) {
#    ok defined $n->pub_date;
#    ok defined $n->create_date;
#    ok defined $n->update_date;
#}
#
## briefly test some Bio::Tree::NodeI methods
#ok my $ancestor = $n->ancestor;
#is $ancestor->scientific_name, 'Homo';
## unless set explicitly, Bio::Taxon doesn't return anything for
## each_Descendent; must ask the database directly
#ok my @children = $ancestor->db_handle->each_Descendent($ancestor);
#cmp_ok @children, '>', 0;
#
#sleep(3) if $db eq $db_entrez;
#
## do some trickier things...
#ok my $n2 = $db->get_Taxonomy_Node('89593');
#is $n2->scientific_name, 'Craniata';
#
## briefly check we can use some Tree methods
#my $tree = Bio::Tree::Tree->new();
#is $tree->get_lca($n, $n2)->scientific_name, 'Craniata';
#
## get lineage_nodes
#my @nodes = $tree->get_nodes;
#is scalar(@nodes), 0;
#my @lineage_nodes;
#@lineage_nodes = $tree->get_lineage_nodes($n->id); # read ID, only works if nodes have been added to tree
#is scalar @lineage_nodes, 0;
#@lineage_nodes = $tree->get_lineage_nodes($n);     # node object always works
#cmp_ok(scalar @lineage_nodes, '>', 20);
#
## get lineage string
#like($tree->get_lineage_string($n), qr/cellular organisms;Eukaryota/);
#like($tree->get_lineage_string($n,'-'), qr/cellular organisms-Eukaryota/);
#like($tree->get_lineage_string($n2), qr/cellular organisms;Eukaryota/);
#
## can we actually form a Tree and use other Tree methods?
#ok $tree = Bio::Tree::Tree->new(-node => $n);
#cmp_ok($tree->number_nodes, '>', 20);
#cmp_ok(scalar($tree->get_nodes), '>', 20);
#is $tree->find_node(-rank => 'genus')->scientific_name, 'Homo';
#
## check that getting the ancestor still works now we have explitly set the
## ancestor by making a Tree
#is $n->ancestor->scientific_name, 'Homo';
#
#sleep(3) if $db eq $db_entrez;
#
#ok $n = $db->get_Taxonomy_Node('1760');
#is $n->scientific_name, 'Actinobacteria';
#
#sleep(3) if $db eq $db_entrez;
#
## entrez isn't as good at searching as flatfile, so we have to special-case
#my @ids = sort $db->get_taxonids('Chloroflexi');
#is scalar @ids, 2;
#is_deeply \@ids, [200795, 32061];
#
#$id = $db->get_taxonids('Chloroflexi (class)');
#$db eq $db_entrez ? is($id, undef) : is($id, 32061);
#
#@ids = $db->get_taxonids('Rhodotorula');
#cmp_ok @ids, '>=' , 8;
#@ids = $db->get_taxonids('Rhodotorula <Microbotryomycetidae>');
#is @ids, 1;
#is $ids[0], 231509;


# get_lca should work on nodes from different databases
#SKIP: {
#    test_skip(-tests => 9, -requires_networking => 1);
#
#    # check that the result is the same as if we are retrieving from the same DB
#    # flatfile
#    $h_flat = $db_flatfile->get_taxon(-name => 'Homo');
#    my $h_flat2 = $db_flatfile->get_taxon(-name => 'Homo sapiens');
#    ok my $tree_functions = Bio::Tree::Tree->new();
#    is $tree_functions->get_lca($h_flat, $h_flat2)->scientific_name, 'Homo', 'get_lca() within flatfile db';
#
#    # entrez
#    my $h_entrez;
#    eval { $h_entrez = $db_entrez->get_taxon(-name => 'Homo sapiens');};
#    skip "Unable to connect to entrez database; no network or server busy?", 7 if $@;
#    my $h_entrez2;
#    eval { $h_entrez2 = $db_entrez->get_taxon(-name => 'Homo');};
#    skip "Unable to connect to entrez database; no network or server busy?", 7 if $@;
#    ok $tree_functions = Bio::Tree::Tree->new();
#    is $tree_functions->get_lca($h_entrez, $h_entrez2)->scientific_name, 'Homo', 'get_lca() within entrez db';
#
#    ok $tree_functions = Bio::Tree::Tree->new();
#    # mixing entrez and flatfile
#    TODO:{
#        local $TODO = 'Mixing databases for get_lca() not working, see bug #3416';
#        is $tree_functions->get_lca($h_flat, $h_entrez)->scientific_name, 'Homo', 'get_lca() mixing flatfile and remote db';
#    }
#    # even though the species taxa for Homo sapiens from list and flat databases
#    # have the same internal id, get_lca won't work because they have different
#    # roots and descendents
#    $h_list = $db_list->get_taxon(-name => 'Homo sapiens');
#    is $h_list->ancestor->internal_id, $h_flat->internal_id;
#    ok ! $tree_functions->get_lca($h_flat, $h_list);
#
#    # but we can form a tree with the flat node then remove all the ranks we're
#    # not interested in and try again
#    $tree = Bio::Tree::Tree->new(-node => $h_flat);
#    $tree->splice(-keep_rank => \@ranks);
#    is $tree->get_lca($h_flat, $h_list)->scientific_name, 'Homo';
#}

# Some tests carried over from flatfile and others that would be nice to pass

## ideas from taxonomy2tree.PLS that let us make nice tree, using
## Bio::Tree::TreeFunctionsI methods; this is a weird and trivial example just
## because our test flatfile database only has the full lineage of one species
#undef $tree;
#for my $name ('Human', 'Hominidae') {
#  my $ncbi_id = $db_flatfile->get_taxonid($name);
#  if ($ncbi_id) {
#    my $node = $db_flatfile->get_taxon(-taxonid => $ncbi_id);
#
#    if ($tree) {
#        ok $tree->merge_lineage($node);
#    }
#    else {
#        ok $tree = Bio::Tree::Tree->new(-node => $node);
#    }
#  }
#}
#is $tree->get_nodes, 30;
#$tree->contract_linear_paths;
#my $ids = join(",", map { $_->id } $tree->get_nodes);
#is $ids, '131567,9606';

#SKIP: {
#    test_skip(-tests => 1, -requires_networking => 1);
#    eval {$db_entrez->get_taxon(10090);};
#    skip "Unable to connect to entrez database; no network or server busy?", 1 if $@;
#
#    my $lca = $db_entrez->get_taxon(314146);
#    my @descs = $db_entrez->get_all_Descendents($lca);
#    cmp_ok @descs, '>=', 17;
#}

END {
    #unlink 'taxonomy.sqlite' if (-e 'taxonomy.sqlite');
}

done_testing();
