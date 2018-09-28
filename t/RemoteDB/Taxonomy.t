# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use Bio::Root::Test;

    test_begin(
        -tests            => 152,
        -requires_modules => [qw(DB_File
                                 XML::Twig )]
    );

    use_ok('Bio::DB::Taxonomy');
    use_ok('Bio::Tree::Tree');
}

my $temp_dir = test_output_dir();

# we're actually testing Bio::Taxon and Bio::DB::Taxonomy::* here, not
# Bio::Taxonomy

ok my $db_flatfile = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
    -namesfile => test_input_file('taxdump', 'names.dmp'),
);
isa_ok $db_flatfile, 'Bio::DB::Taxonomy::flatfile';
isa_ok $db_flatfile, 'Bio::DB::Taxonomy';

# By not specifying a '-directory' argument, index files go to a temporary
# folder ($Bio::Root::IO::TEMPDIR, such as 'C:\Users\USER\AppData\Local\Temp'),
# and are implied to be temporary. So test the ability of flatfile->DESTROY to
# remove the temporary index files at object destruction (this also affects files
# in "test_output_dir()", since the folder is created inside the temporary folder)
no warnings qw(once); # silence 'Name "$Bio::Root::IO::TEMPDIR" used only once'
is $db_flatfile->{index_directory}, $Bio::Root::IO::TEMPDIR, 'removal of temporary index files: no -directory';
$db_flatfile->DESTROY;
ok not -e ($db_flatfile->{index_directory} . '/id2names');
ok not -e ($db_flatfile->{index_directory} . '/names2id');
ok not -e ($db_flatfile->{index_directory} . '/nodes');
ok not -e ($db_flatfile->{index_directory} . '/parents');

# Test removal of temporary index files from test_output_dir folder
# (since test_output_dir() =~ m/^$Bio::Root::IO::TEMPDIR/)
ok $db_flatfile = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => $temp_dir,
    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
    -namesfile => test_input_file('taxdump', 'names.dmp'),
    -force     => 1,
);
is $db_flatfile->{index_directory}, $temp_dir, 'removal of temporary index files: test_output_dir()';
$db_flatfile->DESTROY;
ok not -e ($db_flatfile->{index_directory} . '/id2names');
ok not -e ($db_flatfile->{index_directory} . '/names2id');
ok not -e ($db_flatfile->{index_directory} . '/nodes');
ok not -e ($db_flatfile->{index_directory} . '/parents');

# Generate the object (and the files) again for the remaining tests
ok $db_flatfile = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => $temp_dir,
    -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
    -namesfile => test_input_file('taxdump', 'names.dmp'),
    -force     => 1,
);

{
    my $id;
    my $n;

    is $db_flatfile->get_num_taxa, 189;

    $id = $db_flatfile->get_taxonid('Homo sapiens');

    is $id, 9606;

    # easy test on human, try out the main Taxon methods
    ok $n = $db_flatfile->get_taxon(9606);
    is $n->id, 9606;
    is $n->object_id, $n->id;
    is $n->ncbi_taxid, $n->id;
    is $n->parent_id, 9605;
    is $n->rank, 'species';

    is $n->node_name, 'Homo sapiens';
    is $n->scientific_name, $n->node_name;
    is ${$n->name('scientific')}[0], $n->node_name;

    my %common_names = map { $_ => 1 } $n->common_names;
    cmp_ok keys %common_names, '>=', 3, ref($db_flatfile).": common names";
    ok exists $common_names{human};
    ok exists $common_names{man};

    is $n->division, 'Primates';
    is $n->genetic_code, 1;
    is $n->mitochondrial_genetic_code, 2;

    # briefly test some Bio::Tree::NodeI methods
    ok my $ancestor = $n->ancestor;
    is $ancestor->scientific_name, 'Homo';
    # unless set explicitly, Bio::Taxon doesn't return anything for
    # each_Descendent; must ask the database directly
    ok my @children = $ancestor->db_handle->each_Descendent($ancestor);
    cmp_ok @children, '>', 0;

    # do some trickier things...
    ok my $n2 = $db_flatfile->get_Taxonomy_Node('89593');
    is $n2->scientific_name, 'Craniata';

    # briefly check we can use some Tree methods
    my $tree = Bio::Tree::Tree->new();
    is $tree->get_lca($n, $n2)->scientific_name, 'Craniata';

    # get lineage_nodes
    my @nodes = $tree->get_nodes;
    is scalar(@nodes), 0;
    my @lineage_nodes;
    @lineage_nodes = $tree->get_lineage_nodes($n->id); # read ID, only works if nodes have been added to tree
    is scalar @lineage_nodes, 0;
    @lineage_nodes = $tree->get_lineage_nodes($n);     # node object always works
    cmp_ok(scalar @lineage_nodes, '>', 20);

    # get lineage string
    like($tree->get_lineage_string($n), qr/cellular organisms;Eukaryota/);
    like($tree->get_lineage_string($n,'-'), qr/cellular organisms-Eukaryota/);
    like($tree->get_lineage_string($n2), qr/cellular organisms;Eukaryota/);

    # can we actually form a Tree and use other Tree methods?
    ok $tree = Bio::Tree::Tree->new(-node => $n);
    cmp_ok($tree->number_nodes, '>', 20);
    cmp_ok(scalar($tree->get_nodes), '>', 20);
    is $tree->find_node(-rank => 'genus')->scientific_name, 'Homo';

    # check that getting the ancestor still works now we have explitly set the
    # ancestor by making a Tree
    is $n->ancestor->scientific_name, 'Homo';

    ok $n = $db_flatfile->get_Taxonomy_Node('1760');
    is $n->scientific_name, 'Actinobacteria';

    my @ids = sort $db_flatfile->get_taxonids('Chloroflexi');
    is scalar @ids, 2;
    is_deeply \@ids, [200795, 32061];

    $id = $db_flatfile->get_taxonids('Chloroflexi (class)');
    is($id, 32061);

    @ids = $db_flatfile->get_taxonids('Rhodotorula');
    cmp_ok @ids, '>=' , 1;
    # note the locally cached flatfile is out-of-date, but technically
    # correct for testing purposes
    ok grep { $_ == 266791 } @ids;
    ok grep { $_ == 5533 } @ids;
}


# Test the list database

ok my $db_list = Bio::DB::Taxonomy->new(-source => 'list');
isa_ok $db_list, 'Bio::DB::Taxonomy::list';
isa_ok $db_list, 'Bio::DB::Taxonomy';

my @ranks = qw(superkingdom class genus species);
my @h_lineage = ('Eukaryota', 'Mammalia', 'Homo', 'Homo sapiens');
ok $db_list = Bio::DB::Taxonomy->new(
    -source => 'list',
    -names  => \@h_lineage,
    -ranks  => \@ranks,
);
is $db_list->get_num_taxa, 4;

my @taxa;
ok @taxa = map {$db_list->get_taxon(-name=>$_)} @h_lineage;
is_deeply [map {ref($_)} @taxa], [('Bio::Taxon')x4];
is_deeply [map {$_->rank} @taxa], \@ranks, 'Ranks';

@h_lineage = ('Eukaryota', 'Mammalia', 'Homo', 'Homo erectus');
$db_list->add_lineage(-names => \@h_lineage, -ranks => \@ranks);

ok @taxa = map {$db_list->get_taxon(-name=>$_)} @h_lineage;
is_deeply [map {ref($_)} @taxa], [('Bio::Taxon')x4];
is_deeply [map {$_->rank} @taxa], \@ranks, 'Ranks';

# Make a tree
ok my $tree = $db_list->get_tree('Homo sapiens', 'Homo erectus');
isa_ok $tree, 'Bio::Tree::TreeI';
is $tree->number_nodes, 5;
is $tree->total_branch_length, 4;
ok my $node1 = $tree->find_node( -scientific_name => 'Homo sapiens' );
ok my $node2 = $tree->find_node( -scientific_name => 'Homo erectus' );
is $tree->distance($node1, $node2), 2;

ok my $h_list = $db_list->get_taxon(-name => 'Homo sapiens');
ok my $h_flat = $db_flatfile->get_taxon(-name => 'Homo sapiens');

is $h_list->ancestor->scientific_name, 'Homo';

my @names = $h_list->common_names;
is @names, 0;
$h_list->common_names('woman');
@names = $h_list->common_names;
is @names, 1;
@names = $h_flat->common_names;
is @names, 3;

# you can switch to another database when you need more information, which also
# merges information in the node from the two different dbs
$h_list->db_handle($db_flatfile);
@names = $h_list->common_names;
is @names, 4;

# form a tree with the list lineage first, preventing a subsequent database
# change from giving us all those extra ranks
$h_list->db_handle($db_list);
my $ancestors_ancestor = $h_list->ancestor->ancestor;
is $ancestors_ancestor->scientific_name, 'Mammalia';

$tree = Bio::Tree::Tree->new(-node => $h_list);
$h_list->db_handle($db_flatfile);
$ancestors_ancestor = $h_list->ancestor->ancestor;
is $ancestors_ancestor->scientific_name, 'Mammalia';

# or we can get the flatfile database's idea of the ancestors by removing
# ourselves from the tree
is $h_flat->ancestor->ancestor->scientific_name, 'Homo/Pan/Gorilla group';
$h_list->ancestor(undef);
is $h_list->ancestor->ancestor->scientific_name, 'Homo/Pan/Gorilla group';

# get_lca should work on nodes from different databases
SKIP: {
    test_skip(-tests => 9,
              -requires_modules => ['Bio::DB::Taxonomy::entrez'],
              -requires_networking => 1);
    my $db_entrez = Bio::DB::Taxonomy->new(-source => 'entrez');

    # check that the result is the same as if we are retrieving from the same DB
    # flatfile
    $h_flat = $db_flatfile->get_taxon(-name => 'Homo');
    my $h_flat2 = $db_flatfile->get_taxon(-name => 'Homo sapiens');
    ok my $tree_functions = Bio::Tree::Tree->new();
    is $tree_functions->get_lca($h_flat, $h_flat2)->scientific_name, 'Homo', 'get_lca() within flatfile db';

    # entrez
    my $h_entrez;
    eval { $h_entrez = $db_entrez->get_taxon(-name => 'Homo sapiens');};
    skip "Unable to connect to entrez database; no network or server busy?", 7 if $@;
    my $h_entrez2;
    eval { $h_entrez2 = $db_entrez->get_taxon(-name => 'Homo');};
    skip "Unable to connect to entrez database; no network or server busy?", 7 if $@;
    ok $tree_functions = Bio::Tree::Tree->new();
    is $tree_functions->get_lca($h_entrez, $h_entrez2)->scientific_name, 'Homo', 'get_lca() within entrez db';

    ok $tree_functions = Bio::Tree::Tree->new();
    # mixing entrez and flatfile
    TODO:{
        local $TODO = 'Mixing databases for get_lca() not working, see bug #3416';
        is $tree_functions->get_lca($h_flat, $h_entrez)->scientific_name, 'Homo', 'get_lca() mixing flatfile and remote db';
    }
    # even though the species taxa for Homo sapiens from list and flat databases
    # have the same internal id, get_lca won't work because they have different
    # roots and descendents
    $h_list = $db_list->get_taxon(-name => 'Homo sapiens');
    is $h_list->ancestor->internal_id, $h_flat->internal_id;
    ok ! $tree_functions->get_lca($h_flat, $h_list);

    # but we can form a tree with the flat node then remove all the ranks we're
    # not interested in and try again
    $tree = Bio::Tree::Tree->new(-node => $h_flat);
    $tree->splice(-keep_rank => \@ranks);
    is $tree->get_lca($h_flat, $h_list)->scientific_name, 'Homo';
}

# ideas from taxonomy2tree.PLS that let us make nice tree, using
# Bio::Tree::TreeFunctionsI methods; this is a weird and trivial example just
# because our test flatfile database only has the full lineage of one species
undef $tree;
for my $name ('Human', 'Hominidae') {
  my $ncbi_id = $db_flatfile->get_taxonid($name);
  if ($ncbi_id) {
    my $node = $db_flatfile->get_taxon(-taxonid => $ncbi_id);

    if ($tree) {
        ok $tree->merge_lineage($node);
    }
    else {
        ok $tree = Bio::Tree::Tree->new(-node => $node);
    }
  }
}
is $tree->get_nodes, 30;
$tree->contract_linear_paths;
my $ids = join(",", map { $_->id } $tree->get_nodes);
is $ids, '131567,9606';

# More thorough tests of merge_lineage
ok my $node = $db_list->get_taxon(-name => 'Eukaryota');
$tree = Bio::Tree::Tree->new(-node => $node);
ok $node = $db_list->get_taxon(-name => 'Homo erectus');
ok $tree->merge_lineage($node);
for my $name ('Eukaryota', 'Mammalia', 'Homo', 'Homo erectus') {
   ok $node = $tree->find_node(-scientific_name => $name);
}


# bug 2461
$db_list = Bio::DB::Taxonomy->new(-source => 'list',
                                  -names => [
(split(/,\s+/, "cellular organisms, Eukaryota, Fungi/Metazoa group,
Metazoa, Eumetazoa, Bilateria, Coelomata, Protostomia, Panarthropoda,
Arthropoda, Mandibulata, Pancrustacea, Hexapoda, Insecta, Dicondylia,
Pterygota, Neoptera, Endopterygota, Diptera, Nematocera, Culicimorpha,
Culicoidea, Culicidae, Anophelinae, Anopheles, Anopheles, Angusticorn,
Anopheles, maculipennis group, maculipennis species complex, Anopheles daciae"))]);

my @taxonids = $db_list->get_taxonids('Anopheles');
is @taxonids, 3, 'List context';

my $taxonid = $db_list->get_taxonids('Anopheles');
isa_ok \$taxonid, 'SCALAR', 'Scalar context';
ok exists { map({$_ => undef} @taxonids) }->{$taxonid};

# but we should still be able to merge in an incomplete lineage of a sister
# species and have the 'tree' remain consistent:

# missing 'no rank' Anopheles
$db_list->add_lineage(-names => [
(split(/,\s+/, "Anophelinae, Anopheles, Anopheles, Angusticorn,
maculipennis group, maculipennis species complex, Anopheles labranchiae"))]);
$node = $db_list->get_taxon(-name => 'Anopheles labranchiae');
is $node->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->scientific_name, 'Anophelinae';
is $node->rank, undef;

# missing 'subgenus' Anopheles
$db_list->add_lineage(-names => [
(split(/,\s+/, "Anophelinae, Anopheles, Angusticorn, Anopheles,
maculipennis group, maculipennis species complex, Anopheles maculipennis"))]);
$node = $db_list->get_taxon(-name => 'Anopheles maculipennis');
is $node->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->scientific_name, 'Anophelinae';

# missing 'no rank' Angusticorn
$db_list->add_lineage(-names => [
(split(/,\s+/, "Anophelinae, Anopheles, Anopheles, Anopheles,
maculipennis group, maculipennis species complex, Anopheles melanoon"))]);
$node = $db_list->get_taxon(-name => 'Anopheles melanoon');
is $node->ancestor->ancestor->ancestor->ancestor->scientific_name, 'Angusticorn';

@taxonids = $db_list->get_taxonids('Anopheles');
is scalar @taxonids, 3;

# bug: duplicate topmost taxa
$db_list = Bio::DB::Taxonomy->new( -source => 'list',
                                   -names => ['Bacteria', 'Tenericutes'] );
$db_list->add_lineage(  -names => ['Bacteria'] );
@taxonids = $db_list->get_taxonids('Bacteria');
is scalar @taxonids, 1;

# Disambiguate between taxa with same name using -names
ok $db_list = Bio::DB::Taxonomy->new( -source => 'list' ), 'DB with ambiguous names';
ok $db_list->add_lineage( -names => ['c__Gammaproteobacteria', 'o__Oceanospirillales', 'f__Alteromonadaceae', 'g__Spongiibacter'] );
ok $db_list->add_lineage( -names => ['c__Gammaproteobacteria', 'o__Alteromonadales'  , 'f__Alteromonadaceae', 'g__Alteromonas'  ] );

ok @taxonids = $db_list->get_taxonids('f__Alteromonadaceae');
is scalar @taxonids, 2; # multiple taxa would match using $db_list->get_taxon(-name => 'f__Alteromonadaceae')

ok $node = $db_list->get_taxon( -names => ['c__Gammaproteobacteria', 'o__Alteromonadales'  , 'f__Alteromonadaceae'] );
is $node->ancestor->node_name, 'o__Alteromonadales';
my $iid = $node->internal_id;

ok $node = $db_list->get_taxon( -names => ['c__Gammaproteobacteria', 'o__Oceanospirillales', 'f__Alteromonadaceae'] );
is $node->ancestor->node_name, 'o__Oceanospirillales';
isnt $node->internal_id, $iid;


# More tests with ambiguous names, internal IDs and multiple databases
my ($node3, $node4, $db_list_2);
ok $db_list = Bio::DB::Taxonomy->new( -source => 'list' );
ok $db_list->add_lineage( -names => [ 'o__Enterobacteriales', 'g__Escherichia' ] );
ok $db_list->add_lineage( -names => [ 'o__Pseudomonadales'  , 'g__Pseudomonas' ] );
ok $db_list->add_lineage( -names => [ 'o__Chroococcales'    , 'g__Microcoleus' ] );
ok $node1 = $db_list->get_taxon( -names => [ 'k__Chroococcales', 'g__Microcoleus' ] );

ok $db_list_2 = Bio::DB::Taxonomy->new( -source => 'list' );
ok $db_list_2->add_lineage( -names => [ 'o__Chroococcales', 'g__Microcoleus' ] );
ok $node2 = $db_list_2->get_taxon( -names => [ 'o__Chroococcales', 'g__Microcoleus' ] );

is $node1->scientific_name, 'g__Microcoleus';
is $node2->scientific_name, 'g__Microcoleus'; # same taxon name
isnt $node1->id, $node2->id;                  # but different dbs and hence taxids
is $node1->internal_id, $node1->internal_id;  # but same cross-database internal ID

ok $db_list->add_lineage( -names => [ 'o__Oscillatoriales' , 'g__Microcoleus' ] );
ok $db_list->add_lineage( -names => [ 'o__Acidobacteriales', 'g__Microcoleus' ] );

ok $node1 = $db_list->get_taxon( -names => [ 'o__Chroococcales', 'g__Microcoleus' ] );
ok $node2 = $db_list->get_taxon( -names => [ 'o__Oscillatoriales'  , 'g__Microcoleus' ] );
ok $node3 = $db_list->get_taxon( -names => [ 'o__Acidobacteriales'    , 'g__Microcoleus' ] );
my @nodes = ($node1, $node2, $node3);

is map({$_->id          => undef} @nodes), 6; # 3 distinct taxids
is map({$_->internal_id => undef} @nodes), 6; # 3 distinct iids

ok $db_list->add_lineage( -names => [ 'o__Chroococcales'  , 'g__Microcoleus' ] );
ok $node2 = $db_list->get_taxon( -names => [ 'o__Chroococcales', 'g__Microcoleus' ] );
is $node2->scientific_name, $node1->scientific_name;
is $node2->id, $node1->id;
is $node2->internal_id, $node1->internal_id;
