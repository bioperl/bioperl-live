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


SKIP: {
    test_skip(-tests => 8, -requires_module => 'IO::String');
	
    # test nexus tree parsing
    my $treeio = Bio::TreeIO->new(-format => 'nexus',
				  -verbose => $verbose,
				  -file   => test_input_file('urease.tre.nexus'));
    
    my $tree = $treeio->next_tree;
    ok($tree);
    is($tree->id, 'PAUP_1');
    is($tree->get_leaf_nodes, 6);
    my ($node) = $tree->find_node(-id => 'Spombe');
    is($node->branch_length,0.221404);
    
    # test nexus MrBayes tree parsing
    $treeio = Bio::TreeIO->new(-format => 'nexus',
			       -file   => test_input_file('adh.mb_tree.nexus'));
    
    $tree = $treeio->next_tree;
	my $ct = 1; 
    ok($tree);
    is($tree->id, 'rep.1');
    is($tree->get_leaf_nodes, 54);
    ($node) = $tree->find_node(-id => 'd.madeirensis');
    is($node->branch_length,0.039223);
	while ($tree = $treeio->next_tree) {
		$ct++;
	}
	is($ct,13,'bug 2356');
}


# bug #1854
# process no-newlined tree
my $treeio = Bio::TreeIO->new(-format => 'nexus',
			      -verbose => $verbose,
			      -file   => test_input_file('tree_nonewline.nexus'));

my $tree = $treeio->next_tree;
ok($tree);
ok($tree->find_node('TRXHomo'));

# bug #2205
# process trees with node IDs containing spaces
$treeio = Bio::TreeIO->new(-format => 'nexus',
			   -verbose => $verbose,
			   -file   => test_input_file('spaces.nex'));

$tree = $treeio->next_tree;

my @nodeids = ("'Allium drummondii'", "'Allium cernuum'",'A.cyaneum');

ok($tree);
for my $node ($tree->get_leaf_nodes) {
	is($node->id, shift @nodeids);		
}

# bug #2221
# process tree with names containing quoted commas

$tree = $treeio->next_tree;

@nodeids = ("'Allium drummondii, USA'", "'Allium drummondii, Russia'",'A.cyaneum');

ok($tree);
for my $node ($tree->get_leaf_nodes) {
	is($node->id, shift @nodeids);		
}

# bug #2221
# process tree with names containing quoted commas on one line

$tree = $treeio->next_tree;

@nodeids = ("'Allium drummondii, Russia'", "'Allium drummondii, USA'",'A.cyaneum');

ok($tree);
for my $node ($tree->get_leaf_nodes) {
	is($node->id, shift @nodeids);		
}
