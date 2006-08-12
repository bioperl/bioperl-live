# This is -*-Perl-*- code
# $Id$

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
my $error;

BEGIN {
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	eval {
		require Bio::DB::Taxonomy;
        require Bio::Tree::Tree;
		require XML::Twig;
	};
	if ( $@ ) {
		$error = 1;
		warn "Unable to run tests because XML::Twig is not installed\n";
	}
	$NUMTESTS = 96;
	$error = 0;
	plan tests => $NUMTESTS;
}

END {
    unlink("t/data/taxdump/nodes");
    unlink("t/data/taxdump/parents");
    unlink("t/data/taxdump/id2names");
    unlink("t/data/taxdump/names2id");
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Unable to complete Taxonomy tests',1);
	}
}

if( $error ==  1 ) {
    exit(0);
}

# we're actually testing Bio::Taxon and Bio::DB::Taxonomy::* here, not
# Bio::Taxonomy

my $db_entrez = new Bio::DB::Taxonomy(-source => 'entrez');
ok $db_entrez;

my $dir = "t/data/taxdump";
my $db_flatfile = new Bio::DB::Taxonomy(-source => 'flatfile',
                               -directory => $dir,
                               -nodesfile => Bio::Root::IO->catfile('t','data','taxdump','nodes.dmp'),
                               -namesfile => Bio::Root::IO->catfile('t','data','taxdump','names.dmp'),
                               -force => 1);
ok $db_flatfile;

my $n;
foreach my $db ($db_entrez, $db_flatfile) {
    my $id;
    eval {
        $id = $db->get_taxonid('Homo sapiens');
    };
    if ($@) {
        for (1..32) {
            skip(1, 'Skip unable to connect to entrez database; no network or server busy?');
        }
        next;
    }
    ok $id, 9606;
    
    # easy test on human, try out the main Taxon methods
    ok $n = $db->get_taxon(9606);
    ok $n->id, 9606;
    ok $n->object_id, $n->id;
    ok $n->ncbi_taxid, $n->id;
    ok $n->parent_id, 9605;
    ok $n->rank, 'species';
    
    ok $n->node_name, 'Homo sapiens';
    ok $n->scientific_name, $n->node_name;
    ok ${$n->name('scientific')}[0], $n->node_name;
    
    my %common_names = map { $_ => 1 } $n->common_names;
    ok keys %common_names, 2;
    ok exists $common_names{human};
    ok exists $common_names{man};
    
    ok $n->division, 'Primates';
    ok $n->genetic_code, 1;
    ok $n->mitochondrial_genetic_code, 2;
    # these are entrez-only, data not available in dmp files
    if ($db eq $db_entrez) {
        ok defined $n->pub_date;
        ok defined $n->create_date;
        ok defined $n->update_date;
    }
    
    # briefly test some Bio::Tree::NodeI methods
    ok my $ancestor = $n->ancestor;
    ok $ancestor->scientific_name, 'Homo';
	# unless set explicitly, Bio::Taxon doesn't return anything for
	# each_Descendent; must ask the database directly
    ok my @children = $ancestor->db_handle->each_Descendent($ancestor); 
    ok @children > 0;
    
    sleep(3) if $db eq $db_entrez;
    
    # do some trickier things...
    ok my $n2 = $db->get_Taxonomy_Node('89593');
    ok $n2->scientific_name, 'Craniata';
    
    # briefly check we can use some Tree methods
    my $tree = new Bio::Tree::Tree();
    ok ($tree->get_lca($n, $n2)->scientific_name, 'Craniata');
    
    # can we actually form a Tree and use other Tree methods?
    ok $tree = new Bio::Tree::Tree(-node => $n);
    ok $tree->number_nodes, 30;
    ok $tree->get_nodes, 30;
    ok $tree->find_node(-rank => 'genus')->scientific_name, 'Homo';
    
    # check that getting the ancestor still works now we have explitly set the
    # ancestor by making a Tree
    ok $n->ancestor->scientific_name, 'Homo';
    
    sleep(3) if $db eq $db_entrez;
    
    ok $n = $db->get_Taxonomy_Node('1760');
    ok $n->scientific_name, 'Actinobacteria';
    
    sleep(3) if $db eq $db_entrez;
    
    # entrez isn't as good at searching as flatfile, so we have to special-case
    my @ids = $db->get_taxonids('Chloroflexi');
    $db eq $db_entrez ? (ok @ids, 1) : (ok @ids, 2);
    $id = $db->get_taxonids('Chloroflexi (class)');
    ok $id, 32061;
    
    @ids = $db->get_taxonids('Rhodotorula');
    ok @ids, 8;
    @ids = $db->get_taxonids('Rhodotorula <Microbotryomycetidae>');
    ok @ids, 1;
    ok $ids[0], 231509;
}

# Test the list database
my @ranks = qw(superkingdom class genus species);
my @h_lineage = ('Eukaryota', 'Mammalia', 'Homo', 'Homo sapiens');
my $db_list = new Bio::DB::Taxonomy(-source => 'list', -names => \@h_lineage,
                                                       -ranks => \@ranks);
ok $db_list;

ok my $h_list = $db_list->get_taxon(-name => 'Homo sapiens');
ok my $h_flat = $db_flatfile->get_taxon(-name => 'Homo sapiens');

ok $h_list->ancestor->scientific_name, 'Homo';

my @names = $h_list->common_names;
ok @names, 0;
$h_list->common_names('woman');
@names = $h_list->common_names;
ok @names, 1;
@names = $h_flat->common_names;
ok @names, 2;

# you can switch to another database when you need more information, which also
# merges information in the node from the two different dbs
$h_list->db_handle($db_flatfile);
@names = $h_list->common_names;
ok @names, 3;

# form a tree with the list lineage first, preventing a subsequent database
# change from giving us all those extra ranks
$h_list->db_handle($db_list);
my $ancestors_ancestor = $h_list->ancestor->ancestor;
ok $ancestors_ancestor->scientific_name, 'Mammalia';

my $tree = new Bio::Tree::Tree(-node => $h_list);
$h_list->db_handle($db_flatfile);
$ancestors_ancestor = $h_list->ancestor->ancestor;
ok $ancestors_ancestor->scientific_name, 'Mammalia';

# or we can get the flatfile database's idea of the ancestors by removing
# ourselves from the tree
ok $h_flat->ancestor->ancestor->scientific_name, 'Homo/Pan/Gorilla group';
$h_list->ancestor(undef);
ok $h_list->ancestor->ancestor->scientific_name, 'Homo/Pan/Gorilla group';

# get_lca should work on nodes from different databases
$h_flat = $db_flatfile->get_taxon(-name => 'Homo');
ok my $h_entrez = $db_entrez->get_taxon(-name => 'Homo sapiens');
ok my $tree_functions = new Bio::Tree::Tree();
ok ($tree_functions->get_lca($h_flat, $h_entrez)->scientific_name, 'Homo');

# even though the species taxa for Homo sapiens from list and flat databases
# have the same internal id, get_lca won't work because they have different
# roots and descendents
$h_list = $db_list->get_taxon(-name => 'Homo sapiens');
ok $h_list->ancestor->internal_id, $h_flat->internal_id;
ok ! $tree_functions->get_lca($h_flat, $h_list);

# but we can form a tree with the flat node then remove all the ranks we're
# not interested in and try again
$tree = new Bio::Tree::Tree(-node => $h_flat);
$tree->splice(-keep_rank => \@ranks);
ok ($tree->get_lca($h_flat, $h_list)->scientific_name, 'Homo');

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
		ok $tree = new Bio::Tree::Tree(-node => $node);
    }
  }
}
ok $tree->get_nodes, 30;
$tree->contract_linear_paths;
my $ids = join(",", map { $_->id } $tree->get_nodes);
ok $ids, '131567,9606';