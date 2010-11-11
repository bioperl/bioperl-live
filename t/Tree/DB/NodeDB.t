# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 32);
    
    use_ok('Bio::DB::Tree::Node');
    use_ok('Bio::DB::Tree::Store');
}

my $verbose = test_debug();

# create a db-aware node by hand
my $node = Bio::DB::Tree::Node->new(-node_id => 1);
ok($node,'node created');
isa_ok($node, "Bio::Tree::NodeI");
is($node->node_id, 1,'node id set manually');

# create a simple store
my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
				    -create => 1,
                                    -dsn    => 'dbname=test_tree.idx');
isa_ok($dbh, "Bio::DB::Tree::Store");

# test creating / fetching nodes in the store
my $pk = $dbh->insert_node({'-id' => 'a random label'});
ok($pk,"node created in the database");
$node = $dbh->fetch_node($pk);
ok($node,"found node previously created");
ok(!$node->_dirty,"node is clean");
isa_ok($node,"Bio::Tree::NodeI");
isa_ok($node,"Bio::DB::Tree::Node");
is($node->node_id, $pk, "primary key matches fetch key");
is($node->id, "a random label", "id matches label we provided");
ok(!defined($node->branch_length),"no branch length");
ok(!$node->_dirty,"node is still clean");

# update node properties, then refetch
$node->id("new label");
$node->branch_length(1.5);
ok($node->_dirty, "node is no longer clean");
ok($node->save,"success updating in database");
$node = $dbh->fetch_node($pk);
is($node->id, "new label", "id matches updated label");
is($node->branch_length,1.5,"branch length matches update");

# create a hierarchy of nodes: ((A:2.12,B:3)I1:1,C:6)
$pk = $dbh->insert_node({-id => 'root',-branch_length => 0});
my $root = $dbh->fetch_node($pk);
ok($root,"created and fetched root");
is($root->branch_length,0,'Branch Length for Root');
$pk = $dbh->insert_node({-parent => $root->node_id,
                         -id => 'I1', -branch_length => 1});
my $I1node = $dbh->fetch_node($pk);
is($I1node->branch_length,1,'Branch Length for Internal as integer');
is($I1node->id, "I1", "label matches for Internal");
$pk = $dbh->insert_node({-parent => $I1node->node_id,
                         -id => 'A', -branch_length => 2.12});
my $Anode = $dbh->fetch_node($pk);
is($Anode->branch_length,2.12,'Branch Length for A as float');
is($Anode->id, "A", "label matches for A");
$pk = $dbh->insert_node({-parent => $I1node->node_id,
                         -id => 'B', -branch_length => 3});
my $Bnode = $dbh->fetch_node($pk);
is($Bnode->id, "B", "label matches for B");
$pk = $dbh->insert_node({-parent => $root->node_id,
                         -id => 'C', -branch_length => 6});
my $Cnode = $dbh->fetch_node($pk);
is($Cnode->id, "C", "label matches for C");
is($Cnode->parent_id,$root->node_id,"parent id for C");
is($Bnode->parent_id,$I1node->node_id,"parent id for B");
is($Anode->parent_id,$Bnode->parent_id,"A and B have same parent");

is($Anode->ancestor->node_id,$I1node->node_id,"A's ancestor is I1");
is($Anode->ancestor->node_id,$Bnode->parent_id,"A's ancestor is B's parent");

#my $tree = Bio::DB::Tree::Tree->new(-root_node => $root,
#				    -store     => $dbh);
				    
#my $out = Bio::TreeIO->new(-format => 'newick' );
#$out->write_tree($tree);

done_testing();
unlink('test_tree.idx');
exit;
