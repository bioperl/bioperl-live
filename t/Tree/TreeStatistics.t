# -*-Perl-*- Test Harness script for Bioperl
# $Id: RandomTreeFactory.t 11525 2007-06-27 10:16:38Z sendu $

use strict;
use FindBin qw/$RealBin/;
use lib "$RealBin/../../lib";

BEGIN { 
    use Bio::Root::Test;
    
    test_begin(-tests => 44);

    use_ok('Bio::TreeIO');
    use_ok('Bio::Tree::Statistics');
}


use Data::Dumper;
my $in = Bio::TreeIO->new(-format => 'nexus',
                          -file   => test_input_file('traittree.nexus'));
my $tree = $in->next_tree;
my $node = $tree->find_node(-id => 'N14');

# Create some "bootstrap" trees for the next couple of tests
my @bs_trees = (1) x 10; # earmark the memory but clone the tree in the next loop
# Alter the the trees so that they end up with less
# than 100% support

for(my $bsTreeIndex=0; $bsTreeIndex < @bs_trees; $bsTreeIndex+=1){

  $bs_trees[$bsTreeIndex] = $tree->clone;
  my @bsLeaf = sort {$a->id <=> $b->id } grep{$_->is_Leaf} $bs_trees[$bsTreeIndex]->get_nodes;
  # Mix the first node with the $bsTreeIndex node
  my $leafIndex = $bsTreeIndex % int(scalar(@bsLeaf)/2); # only messing with 1/2 the leaves
  my($name1,$name2) = ($bsLeaf[0]->id, $bsLeaf[$leafIndex]->id);
  $bsLeaf[0]         ->id($name2);
  $bsLeaf[$leafIndex]->id($name1);
  
  # Mess with a second taxon
  my $leafIndex = $bsTreeIndex % scalar(@bsLeaf); # mess with all leaves
  my($name3,$name4) = ($bsLeaf[-1]->id, $bsLeaf[$leafIndex]->id);
  $bsLeaf[-1]        ->id($name4);
  $bsLeaf[$leafIndex]->id($name3);
}

my $stats = Bio::Tree::Statistics->new();
is $stats->cherries($tree), 8, 'cherries';
is $stats->cherries($tree, $node), 4, 'cherries';

subtest 'transfer-bootstrap-expectation (experimental)' => sub{
  plan tests=>15;
  my %expectation = (''=>100, N1=>27, N2=>82, N3=>64, N4=>82, N5=>82, N8=>82, N6=>82, N7=>82, N9=>100, N10=>91, N11=>100, N12=>9, N13=>55, N14=>82);
  my $bs_tree  = $stats->transfer_bootstrap_expectation(\@bs_trees, $tree);
  my @node = $bs_tree->get_nodes;
  for(my $i=0;$i<@node;$i++){
    next if($node[$i]->is_Leaf);
    is $node[$i]->bootstrap , $expectation{$node[$i]->id}, "Testing TBE for node ".$node[$i]->id;
  }
};


subtest 'assess_bootstrap' => sub{
  plan tests=>15;
  my %expectation = (''=>100, N1=>20, N2=>80, N3=>20, N4=>80, N5=>80, N8=>80, N6=>60, N7=>20, N9=>100, N10=>80, N11=>100, N12=>9, N13=>55, N14=>20);
  my $bs_tree  = $stats->assess_bootstrap(\@bs_trees, $tree);
  my @node = $bs_tree->get_nodes;
  for(my $i=0;$i<@node;$i++){
    next if($node[$i]->is_Leaf);
    is $node[$i]->bootstrap, $expectation{$node[$i]->id}, "Testing bootstrap for node ".$node[$i]->id;
  }
};

# traits
my $key = $tree->add_trait(test_input_file('traits.tab'), 4);
is $key, undef, 'read traits'; # exceeded column number

$key = $tree->add_trait(test_input_file('traits.tab'), 2, 1);
is $key, 'disp', "Add traits in second column and ignore missing"; # one leaf has a missing trait value, but ignore it

$key = $tree->add_trait(test_input_file('traits.tab'), 3);
is $key, 'intermediate', "Add traits in third column";

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

