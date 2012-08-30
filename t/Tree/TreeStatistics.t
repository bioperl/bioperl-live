# -*-Perl-*- Test Harness script for Bioperl
# $Id: RandomTreeFactory.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 42);

    use_ok('Bio::TreeIO');
    use_ok('Bio::Tree::Statistics');
}


use Data::Dumper;
my $in = Bio::TreeIO->new(-format => 'nexus',
                          -file   => test_input_file('traittree.nexus'));
my $tree = $in->next_tree;
my $node = $tree->find_node(-id => 'N14');


my $stats = Bio::Tree::Statistics->new();
is $stats->cherries($tree), 8, 'cherries';
is $stats->cherries($tree, $node), 4, 'cherries';

# traits
my $key = $tree->add_trait(test_input_file('traits.tab'), 4);
is $key, undef, 'read traits'; # exceeded column number

$key = $tree->add_trait(test_input_file('traits.tab'), 2, 1);
is $key, 'disp'; # one leaf has a missing trait value, but ignore it

$key = $tree->add_trait(test_input_file('traits.tab'), 3);
is $key, 'intermediate';

is $stats->ps($tree, $key), 4, 'parsimony score';
is $stats->ps($tree, $key, $node), 1, 'subtree parsimony score';

my $node_i = $tree->find_node(-id => 'N10');
my @values = sort $node_i->get_tag_values('ps_trait');
ok eq_set (\@values, ['red', 'blue']), 'ps value';

is $stats->fitch_down($tree), 1, 'fitch_down';
is $node_i->get_tag_values('ps_trait'), 'red', 'ps value after fitch_down';



$node_i = $tree->find_node(-id => '2'); # leaf
is $stats->persistence($tree, $node_i), 1, 'persistence of a leaf';

$node_i = $tree->find_node(-id => 'N1');
is $stats->persistence($tree, $node_i), 1, 'persistence of an internal node value ';

$node_i = $tree->find_node(-id => 'N13');
is $stats->persistence($tree, $node_i), 3,  'persistence of an internal node value';

$node_i = $tree->find_node(-id => 'N6');
is $stats->persistence($tree, $node_i), 2,  'persistence of an internal node value';

my $value;

$node_i = $tree->find_node(-id => '1');
is $stats->count_subclusters($tree, $node_i), 0,  'leaf node: number of clusters = 0 ';

$node_i = $tree->find_node(-id => 'N3');
is $stats->count_subclusters($tree, $node_i), 1,  'number of clusters ';

$node_i = $tree->find_node(-id => 'N14');
is $stats->count_subclusters($tree, $node_i), 1,  'number of clusters ';

$node_i = $tree->find_node(-id => 'N7');
is $stats->count_subclusters($tree, $node_i), 2,  'number of clusters ';



$node_i = $tree->find_node(-id => 'N12');
is $stats->count_leaves($tree, $node_i), 2,  'number of leaves in phylotype ';

$node_i = $tree->find_node(-id => 'N13');
is $stats->count_leaves($tree, $node_i), 4,  'number of leaves in phylotype ';

$node_i = $tree->find_node(-id => 'N14');
is $stats->count_leaves($tree, $node_i), 6,  'number of leaves in phylotype ';

$node_i = $tree->find_node(-id => 'N7');
is $stats->count_leaves($tree, $node_i), 6,  'number of leaves in phylotype ';



$node_i = $tree->find_node(-id => 'N4');
is $stats->phylotype_length($tree, $node_i), 1,  'phylotype length';

$node_i = $tree->find_node(-id => 'N6');
is $stats->phylotype_length($tree, $node_i), 5,  'phylotype length';

$node_i = $tree->find_node(-id => 'N7');
is $stats->phylotype_length($tree, $node_i), 12,  'phylotype length';

$node_i = $tree->find_node(-id => 'N13');
is $stats->phylotype_length($tree, $node_i), 6,  'phylotype length';

$node_i = $tree->find_node(-id => 'N14');
is $stats->phylotype_length($tree, $node_i), 11,  'phylotype length';


$node_i = $tree->find_node(-id => 'N4');
is $stats->sum_of_leaf_distances($tree, $node_i), 1,  'sum of leaf distances';

$node_i = $tree->find_node(-id => 'N6');
is $stats->sum_of_leaf_distances($tree, $node_i), 6,  'sum of leaf distances';

$node_i = $tree->find_node(-id => 'N7');
is $stats->sum_of_leaf_distances($tree, $node_i), 18,  'sum of leaf distances';

$node_i = $tree->find_node(-id => 'N13');
is $stats->sum_of_leaf_distances($tree, $node_i), 8,  'sum of leaf distances';

$node_i = $tree->find_node(-id => 'N14');
is $stats->sum_of_leaf_distances($tree, $node_i), 18,  'sum of leaf distances';



is sprintf ("%.3f", $stats->genetic_diversity($tree, $node_i)), '3.000',  'genetic diversity'; 

is sprintf ("%.3f", $stats->statratio($tree, $node_i)), '0.333',  'separation'; 


is $stats->ai($tree, $key), 0.628906, 'association index';
is $stats->ai($tree, $key, $node), 0.062500, 'subtree association index';

my $mc = $stats->mc($tree, $key);
is ($mc->{blue}, 2, 'monophyletic clade size');
is ($mc->{red}, 4, 'monophyletic clade size');
$node = $tree->find_node(-id => 'N10');
$mc = $stats->mc($tree, $key, $node);
is ($mc->{blue}, 2, 'monophyletic clade size');
is ($mc->{red}, 2, 'monophyletic clade size');


