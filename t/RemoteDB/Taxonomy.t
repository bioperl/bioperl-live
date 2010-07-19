# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 103,
			   -requires_module => 'XML::Twig');
	
	use_ok('Bio::DB::Taxonomy');
	use_ok('Bio::Tree::Tree');
}

my $temp_dir = test_output_dir();

# we're actually testing Bio::Taxon and Bio::DB::Taxonomy::* here, not
# Bio::Taxonomy

ok my $db_entrez = Bio::DB::Taxonomy->new(-source => 'entrez');

ok my $db_flatfile = Bio::DB::Taxonomy->new(-source => 'flatfile',
                               -directory => $temp_dir,
                               -nodesfile => test_input_file('taxdump', 'nodes.dmp'),
                               -namesfile => test_input_file('taxdump','names.dmp'),
                               -force => 1);

my $n;
foreach my $db ($db_entrez, $db_flatfile) {
    SKIP: {
		test_skip(-tests => 38, -requires_networking => 1) if $db eq $db_entrez;
        my $id;
        eval { $id = $db->get_taxonid('Homo sapiens');};
        skip "Unable to connect to entrez database; no network or server busy?", 38 if $@;
        
        is $id, 9606;
        
        # easy test on human, try out the main Taxon methods
        ok $n = $db->get_taxon(9606);
        is $n->id, 9606;
        is $n->object_id, $n->id;
        is $n->ncbi_taxid, $n->id;
        is $n->parent_id, 9605;
        is $n->rank, 'species';
        
        is $n->node_name, 'Homo sapiens';
        is $n->scientific_name, $n->node_name;
        is ${$n->name('scientific')}[0], $n->node_name;
        
        my %common_names = map { $_ => 1 } $n->common_names;
        is keys %common_names, 3, ref($db).": common names";
        ok exists $common_names{human};
        ok exists $common_names{man};
        
        is $n->division, 'Primates';
        is $n->genetic_code, 1;
        is $n->mitochondrial_genetic_code, 2;
        # these are entrez-only, data not available in dmp files
        if ($db eq $db_entrez) {
            ok defined $n->pub_date;
            ok defined $n->create_date;
            ok defined $n->update_date;
        }
        
        # briefly test some Bio::Tree::NodeI methods
        ok my $ancestor = $n->ancestor;
        is $ancestor->scientific_name, 'Homo';
        # unless set explicitly, Bio::Taxon doesn't return anything for
        # each_Descendent; must ask the database directly
        ok my @children = $ancestor->db_handle->each_Descendent($ancestor); 
        ok @children > 0;
        
        sleep(3) if $db eq $db_entrez;
        
        # do some trickier things...
        ok my $n2 = $db->get_Taxonomy_Node('89593');
        is $n2->scientific_name, 'Craniata';
        
        # briefly check we can use some Tree methods
        my $tree = Bio::Tree::Tree->new();
        is $tree->get_lca($n, $n2)->scientific_name, 'Craniata';
        
        # can we actually form a Tree and use other Tree methods?
        ok $tree = Bio::Tree::Tree->new(-node => $n);
        is $tree->number_nodes, 30;
        is $tree->get_nodes, 30;
        is $tree->find_node(-rank => 'genus')->scientific_name, 'Homo';
        
        # check that getting the ancestor still works now we have explitly set the
        # ancestor by making a Tree
        is $n->ancestor->scientific_name, 'Homo';
        
        sleep(3) if $db eq $db_entrez;
        
        ok $n = $db->get_Taxonomy_Node('1760');
        is $n->scientific_name, 'Actinobacteria';
        
        sleep(3) if $db eq $db_entrez;
        
        # entrez isn't as good at searching as flatfile, so we have to special-case
        my @ids = $db->get_taxonids('Chloroflexi');
        $db eq $db_entrez ? (is @ids, 1) : (is @ids, 2);
        $id = $db->get_taxonids('Chloroflexi (class)');
        is $id, 32061;
        
        @ids = $db->get_taxonids('Rhodotorula');
        cmp_ok @ids, '>=' , 8;
        @ids = $db->get_taxonids('Rhodotorula <Microbotryomycetidae>');
        is @ids, 1;
        is $ids[0], 231509;
    }
}

# Test the list database
my @ranks = qw(superkingdom class genus species);
my @h_lineage = ('Eukaryota', 'Mammalia', 'Homo', 'Homo sapiens');
my $db_list = Bio::DB::Taxonomy->new(-source => 'list', -names => \@h_lineage,
                                                        -ranks => \@ranks);
ok $db_list;

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

my $tree = Bio::Tree::Tree->new(-node => $h_list);
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
    test_skip(-tests => 5, -requires_networking => 1);
    $h_flat = $db_flatfile->get_taxon(-name => 'Homo');
    my $h_entrez;
    eval { $h_entrez = $db_entrez->get_taxon(-name => 'Homo sapiens');};
    skip "Unable to connect to entrez database; no network or server busy?", 5 if $@;
    
    ok my $tree_functions = Bio::Tree::Tree->new();
    is $tree_functions->get_lca($h_flat, $h_entrez)->scientific_name, 'Homo';
    
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
		$tree->merge_lineage($node);
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

# we can recursively fetch all descendents of a taxon
SKIP: {
    test_skip(-tests => 1, -requires_networking => 1);
    eval {$db_entrez->get_taxon(10090);};
    skip "Unable to connect to entrez database; no network or server busy?", 1 if $@;
    
    my $lca = $db_entrez->get_taxon(314146);
    my @descs = $db_entrez->get_all_Descendents($lca);
    cmp_ok @descs, '>=', 17;
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
is @taxonids, 3;

# but we should still be able to merge in an incomplete lineage of a sister
# species and have the 'tree' remain consistent:

# missing 'no rank' Anopheles
$db_list->add_lineage(-names => [
(split(/,\s+/, "Anophelinae, Anopheles, Anopheles, Angusticorn,
maculipennis group, maculipennis species complex, Anopheles labranchiae"))]);
my $node = $db_list->get_taxon(-name => 'Anopheles labranchiae');
is $node->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->ancestor->scientific_name, 'Anophelinae';

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
is @taxonids, 3;
