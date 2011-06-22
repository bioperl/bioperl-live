#-*-perl-*-
# $Id$

use strict;
use Bio::Root::Test;
test_begin( -tests=>46,
	    -requires_modules => [qw(Bio::Phylo)]);

use_ok( 'Bio::Tree::Tree' );
use_ok( 'Bio::TreeIO' );
use_ok('Bio::TreeIO::nexml'); # checks that your module is there and loads ok

diag("WARNING: NeXML parsing for NeXML v0.9 is currently very experimental support");
#Read in Data
ok( my $TreeStream = Bio::TreeIO->new(-file => test_input_file('nexml','trees.nexml.xml'), -format => 'nexml') );

	#checking first tree object
	ok( my $tree_obj = $TreeStream->next_tree(), 'tree obj read' );
	isa_ok($tree_obj, 'Bio::Tree::Tree');
	is( $tree_obj->get_root_node()->id(), 'n1', "root node");
	my @nodes = $tree_obj->get_nodes();
	is( @nodes, 9, "number of nodes");
	ok ( my $node7 = $tree_obj->find_node('n7') );
	is( $node7->branch_length, 0.3247, "branch length");
	is( $node7->ancestor->id, 'n3', 'anc id');
	is( $node7->ancestor->branch_length, '0.34534', 'anc bl');
	#Check leaf nodes and taxa
	my %expected_leaves = (
							'n8'	=>	'bird',
							'n9'	=>	'worm',
							'n5'	=>	'dog',
							'n6'	=>	'mouse',
							'n2'	=>	'human'
	);
	ok( my @leaves = $tree_obj->get_leaf_nodes() );
	is( @leaves, 5, "number of leaf nodes");
	foreach my $leaf (@leaves) {
		my $leafID = $leaf->id();
		ok( exists $expected_leaves{$leaf->id()}, "$leafID exists"  );
		is( $leaf->get_tag_values('taxon'), $expected_leaves{$leaf->id()}, "$leafID taxon");
	}
	
	
#Write data
diag('Begin tests for writing tree files');
my $outdata = test_output_file();
ok( my $outTreeStream = Bio::TreeIO->new(-file => ">$outdata", -format => 'nexml'), 'out stream');
ok( $outTreeStream->write_tree($tree_obj), 'write tree');
close($outdata);

#Read in the out file to test roundtrip
my $inTreeStream = Bio::TreeIO->new(-file => $outdata, -format => 'nexml');
	
	#checking first tree object
	ok($tree_obj = $inTreeStream->next_tree(), 'read tree obj (rt)' );
	isa_ok($tree_obj, 'Bio::Tree::Tree');
	is( $tree_obj->get_root_node()->id(), 'n1', "root node (rt)");
	my @outnodes = $tree_obj->get_nodes();
	is( @outnodes, 9, "number of nodes (rt)");
	ok ( $node7 = $tree_obj->find_node('n7') );
	is( $node7->branch_length, 0.3247, "branch length (rt)");
	is( $node7->ancestor->id, 'n3','anc id (rt)');
	is( $node7->ancestor->branch_length, '0.34534', 'anc bl (rt)');
	
	#Check leaf nodes and taxa
	ok( my @outleaves = $tree_obj->get_leaf_nodes() );
	is( @outleaves, 5, "number of leaf nodes (rt)");
	foreach my $leaf (@outleaves)
	{
		my $leafID = $leaf->id();
		ok( exists $expected_leaves{$leaf->id()}, "$leafID exists (rt)"  );
		is( $leaf->get_tag_values('taxon'), $expected_leaves{$leaf->id()}, "$leafID taxon (rt)");
	}

	
