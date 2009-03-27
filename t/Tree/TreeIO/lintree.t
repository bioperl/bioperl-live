# -*-Perl-*- Test Harness script for Bioperl
# $Id: TreeIO.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 24);
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

# try lintree parsing
my $treeio = Bio::TreeIO->new(-format => 'lintree',
			      -file   => test_input_file('crab.njb'));

my ($tree, @leaves, $node);
while( $tree = $treeio->next_tree ) {

    isa_ok($tree, 'Bio::Tree::TreeI');

    my @nodes = $tree->get_nodes;

    @leaves = $tree->get_leaf_nodes;
    is(@leaves, 13);
#\maj    is(@nodes, 25);
    is(@nodes, 24);
    ($node) = $tree->find_node(-id => '18');
    ok($node);
    is($node->id, '18');
    is($node->branch_length, '0.030579');
    is($node->bootstrap, 998);
}

$treeio = Bio::TreeIO->new(-format => 'lintree',
			   -file   => test_input_file('crab.nj'));

$tree = $treeio->next_tree;

isa_ok($tree, 'Bio::Tree::TreeI');

my @nodes = $tree->get_nodes;
@leaves = $tree->get_leaf_nodes;
is(@leaves, 13);
#/maj is(@nodes, 25);
is(@nodes, 24); #/maj
($node) = $tree->find_node('18');
is($node->id, '18');
is($node->branch_length, '0.028117');

($node) = $tree->find_node(-id => 'C-vittat');
is($node->id, 'C-vittat');
is($node->branch_length, '0.087619');
is($node->ancestor->id, '14');

$treeio = Bio::TreeIO->new(-format => 'lintree',
			  -file   => test_input_file('crab.dat.cn'));

$tree = $treeio->next_tree;

isa_ok($tree, 'Bio::Tree::TreeI');

@nodes = $tree->get_nodes;
@leaves = $tree->get_leaf_nodes;
is(@leaves, 13, "Leaf nodes");

#/maj is(@nodes, 25, "All nodes");
is(@nodes, 24, "All nodes"); #/maj
($node) = $tree->find_node('18');
is($node->id, '18');

is($node->branch_length, '0.029044');

($node) = $tree->find_node(-id => 'C-vittat');
is($node->id, 'C-vittat');
is($node->branch_length, '0.097855');
is($node->ancestor->id, '14');
