# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 35);
    
    use_ok('Bio::DB::Tree::Tree');
    use_ok('Bio::DB::Tree::Node');
    use_ok('Bio::DB::Tree::Store');
}

my $verbose = test_debug();

# create a db-aware tree by hand
my $tree = Bio::DB::Tree::Tree->new(-tree_id => 1);
ok($tree,'tree created');
isa_ok($tree, "Bio::Tree::TreeI");
is($tree->tree_id, 1,'tree id set manually');

# create a simple store
my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
				    -create => 1,
                                    -dsn    => 'dbname=test_tree.idx');
isa_ok($dbh, "Bio::DB::Tree::Store");

# test creating / fetching trees in the store
my $pk = $dbh->insert_tree({'-id' => 'tree1'});
ok($pk,"tree created in the database");
$tree = $dbh->fetch_tree($pk);
ok($tree,"found tree previously created");
ok(!$tree->_dirty,"tree is clean");
isa_ok($tree,"Bio::Tree::TreeI");
isa_ok($tree,"Bio::DB::Tree::Tree");
is($tree->tree_id, $pk, "primary key matches fetch key");
is($tree->id, "tree1", "id matches label we provided");
ok($tree->is_rooted, "by default tree is rooted");
ok(! defined($tree->root_node), "no root node");
ok(!$tree->_dirty,"tree is still clean");

# repeat with a root node
my $root = $dbh->insert_node({-id => "tree2 root"});
my $pk2 = $dbh->insert_tree({-id => 'tree2', -root => $root});
ok($pk2,"tree2 created in the database");
my $tree2 = $dbh->fetch_tree($pk2);
ok($tree2,"found tree2 by key");
is($tree2->id, "tree2", "id matches label we provided");
ok($tree2->is_rooted, "by default tree is rooted");
isa_ok($tree2->root_node,"Bio::DB::Tree::Node");
is($tree2->root_node->node_id,$root,"root node is what we created");
is($tree2->root_node->id,"tree2 root","root label is what we created");
ok(!$tree2->_dirty,"tree is still clean");

# update properties, then refetch
$tree->rooted(0);
$tree->root_node($dbh->fetch_node($root));
$tree->root_node->id("mytree root");
$tree->id("mytree");
ok($tree->_dirty, "tree is no longer clean");
ok($tree->root_node->save,"success updating root node in the database");
ok($tree->save,"success updating tree in database");
$tree = $dbh->fetch_tree($tree->tree_id);
is($tree->id, "mytree", "id matches updated label");
ok(!$tree->is_rooted, "not rooted after update");
ok(!$tree->rooted, "not rooted after update");
ok(defined($tree->root_node), "has \"root\" node");
is($tree->root_node->node_id,$root,"root node is what we set it to");
is($tree->root_node->id,"mytree root","root node has updated state");

# create a hierarchy of nodes: ((A:2.12,B:3)I1:1,C:6)
$pk = $dbh->insert_node({-id => 'root',-branch_length => 0});
$root = $dbh->fetch_node($pk);
$pk = $dbh->insert_node({-parent => $root->node_id,
                         -id => 'I1', -branch_length => 1});
my $I1node = $dbh->fetch_node($pk);
$pk = $dbh->insert_node({-parent => $I1node->node_id,
                         -id => 'A', -branch_length => 2.12});
my $Anode = $dbh->fetch_node($pk);
$pk = $dbh->insert_node({-parent => $I1node->node_id,
                         -id => 'B', -branch_length => 3});
my $Bnode = $dbh->fetch_node($pk);
$pk = $dbh->insert_node({-parent => $root->node_id,
                         -id => 'C', -branch_length => 6});
my $Cnode = $dbh->fetch_node($pk);
$tree = $dbh->fetch_tree($dbh->insert_tree({-id => 'test tree', 
                                            -root => $root->node_id}));
ok($tree,"created test tree");
#my $out = Bio::TreeIO->new(-format => 'newick' );
#$out->write_tree($tree);

done_testing();
unlink('test_tree.idx');
exit;
