# -*-Perl-*- Test Harness script for Bioperl
# $Id: TreeIO.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 7);
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

my $treeio = Bio::TreeIO->new(-verbose => $verbose,
			  -file   => test_input_file('test.nhx'));
my $tree;



ok($treeio);
$tree = $treeio->next_tree;

isa_ok($tree, 'Bio::Tree::TreeI');

my @nodes = $tree->get_nodes;
is(@nodes, 13, "Total Nodes");

my $adhy = $tree->find_node('ADHY');
is($adhy->branch_length, 0.1);
is(($adhy->get_tag_values('S'))[0], 'nematode');
is(($adhy->get_tag_values('E'))[0], '1.1.1.1');
