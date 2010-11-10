# -*-Perl-*- Test Harness script for Bioperl

use strict;

use Bio::Root::Test;
use Bio::DB::Tree::Node;
use Bio::DB::Tree::Tree;
use Bio::DB::Tree::Store;
use Bio::TreeIO;

my $node = Bio::DB::Tree::Node->new(-node_id => 1);
ok($node,'node created');
is($node->node_id, 1,'node id set manually');

my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'SQLite',
				    -create => 1,
				    -dsn    => 'dbname=test_tree.idx');

my $rootnode = $dbh->_create_node(undef,'root',0,undef);
ok($rootnode,'root node created');
is($rootnode->branch_length,0,'Branch Length for Root');
my $intnode = $dbh->_create_node($rootnode->node_id,'I1',1,undef);
is($intnode->branch_length,1,'Branch Length for Internal');
my $Anode = $dbh->_create_node($intnode->node_id,'A',2.12,undef);
is($Anode->branch_length,2.12,'Branch Length for A with real');
my $Bnode = $dbh->_create_node($intnode->node_id,'B',3,undef);
my $Cnode = $dbh->_create_node($rootnode->node_id,'A',6,undef);
is($Anode->ancestor->node_id,$intnode->node_id);
is($Anode->ancestor->node_id,$Bnode->parent_id);

my $tree = Bio::DB::Tree::Tree->new(-root_node => $rootnode,
				    -store     => $dbh);
				    
#my $out = Bio::TreeIO->new(-format => 'newick' );
#$out->write_tree($tree);

done_testing();
unlink('test_tree.idx');
exit;
