# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 42);
	
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

my $treeio = Bio::TreeIO->new(-verbose => $verbose,
			     -format => 'nhx',
			     -file   => test_input_file('test.nhx'));
my $tree = $treeio->next_tree;

my @nodes = $tree->find_node('ADH2');
is(@nodes, 2,'Number of nodes that have ADH2 as name');

if( $verbose ) {
    $treeio = Bio::TreeIO->new(-verbose => $verbose,
			      -format => 'nhx',
			      );
    $treeio->write_tree($tree);
    print "nodes are: \n",
    join(", ", map {  $_->id . ":". (defined $_->branch_length ? 
				     $_->branch_length : '' ) } @nodes), "\n";
}

$treeio = Bio::TreeIO->new(-format => 'newick',
			  -file   => test_input_file('test.nh'));
$tree = $treeio->next_tree;


if( $verbose ) { 
    my $out = Bio::TreeIO->new(-format => 'tabtree');
    
    $out->write_tree($tree);
}

my @hADH = ( $tree->find_node('hADH1'),
	     $tree->find_node('hADH2') );
my ($n4) = $tree->find_node('yADH4');

is($tree->is_monophyletic(-nodes    => \@hADH,
			  -outgroup => $n4),1,'Test Monophyly');

my @mixgroup = ( $tree->find_node('hADH1'),
		 $tree->find_node('yADH2'),
		 $tree->find_node('yADH3'),
		 );

my ($iADHX) = $tree->find_node('iADHX');

# test height
is($iADHX->height, 0,'Height');
is($iADHX->depth,0.22,'Depth');
isnt( $tree->is_monophyletic(-nodes   => \@mixgroup,
			    -outgroup=> $iADHX),1, 'non-monophyletic group');

my $in = Bio::TreeIO->new(-format => 'newick',
			 -fh     => \*DATA);
$tree = $in->next_tree;
my ($a,$b,$c,$d) = ( $tree->find_node('A'),
		     $tree->find_node('B'),
		     $tree->find_node('C'),
		     $tree->find_node('D'));

is($tree->is_monophyletic(-nodes => [$b,$c],
			  -outgroup => $d),1, 'B,C are Monophyletic');

is($tree->is_monophyletic(-nodes => [$b,$a],
			  -outgroup => $d),1,'A,B are Monophyletic');

$tree = $in->next_tree;
my ($e,$f,$i);
($a,$b,$c,$d,$e,$f,$i) = ( $tree->find_node('A'),
			   $tree->find_node('B'),
			   $tree->find_node('C'),
			   $tree->find_node('D'),
			   $tree->find_node('E'),
			   $tree->find_node('F'),
			   $tree->find_node('I'),
			   );
isnt( $tree->is_monophyletic(-nodes => [$b,$f],
			    -outgroup => $d),1,'B,F are not Monophyletic' );

is($tree->is_monophyletic(-nodes => [$b,$a],
			  -outgroup => $f),1, 'A,B are Monophyletic');

# test for paraphyly

isnt(  $tree->is_paraphyletic(-nodes => [$a,$b,$c],
			   -outgroup => $d), 1,'A,B,C are not Monophyletic w D as outgroup');

is(  $tree->is_paraphyletic(-nodes => [$a,$f,$e],
			   -outgroup => $i), 1, 'A,F,E are monophyletic with I as outgroup');


# test for rerooting the tree
my $out = Bio::TreeIO->new(-format => 'newick', 
			   -fh => \*STDERR, 
			   -noclose => 1);
$tree = $in->next_tree;
$tree->verbose( -1 ) unless $verbose;
my $node_cnt_orig = scalar($tree->get_nodes);
# reroot on an internal node: should work fine
$a = $tree->find_node('A');
# removing node_count checks because re-rooting can change the
# number of internal nodes (if it is done correctly)
my $total_length_orig = $tree->total_branch_length;
$out->write_tree($tree) if $verbose;
is($tree->reroot($a),1, 'Can re-root with A as outgroup');
$out->write_tree($tree) if $verbose;
is($node_cnt_orig, scalar($tree->get_nodes), 'Count the number of nodes');
my $total_length_new = $tree->total_branch_length;
my $eps = 0.001 * $total_length_new;	# tolerance for checking length
warn("orig total len ", $total_length_orig, "\n") if $verbose;
warn("new  total len ", $tree->total_branch_length,"\n") if $verbose;
# according to retree in phylip these branch lengths actually get larger
# go figure...
#ok(($total_length_orig >= $tree->total_branch_length - $eps)
#   and ($total_length_orig <= $tree->total_branch_length + $eps));
is($tree->get_root_node, $a->ancestor, "Root node is A's ancestor");

# try to reroot on an internal, will result in there being 1 less node
$a = $tree->find_node('C')->ancestor;
$out->write_tree($tree) if $verbose;
is($tree->reroot($a),1, "Can reroot with C's ancsestor");
$out->write_tree($tree) if $verbose;
is($node_cnt_orig, scalar($tree->get_nodes), 'Check to see that node count is correct after an internal node was removed after this re-rooting');
warn("orig total len ", $total_length_orig, "\n") if $verbose;
warn("new  total len ", $tree->total_branch_length,"\n") if $verbose;
cmp_ok($total_length_orig, '>=', $tree->total_branch_length - $eps, 
       'Total original branch length is what it is supposed to be');
cmp_ok($total_length_orig, '<=',$tree->total_branch_length + $eps, 
       'Updated total branch length after the reroot');
is($tree->get_root_node, $a->ancestor, 'Make sure root is really what we asked for');

# try to reroot on existing root: should fail
$a = $tree->get_root_node;
isnt( $tree->reroot($a),1, 'Testing for failed re-rerooting');

# try a more realistic tree
$tree = $in->next_tree;
$a = $tree->find_node('VV');
$node_cnt_orig = scalar($tree->get_nodes);
$total_length_orig = $tree->total_branch_length;
$out->write_tree($tree) if $verbose;
is($tree->reroot($a->ancestor),1, 'Test that rooting succeeded');
$out->write_tree($tree) if $verbose;
is($node_cnt_orig+1, scalar($tree->get_nodes), 'Test that re-rooted tree has proper number of nodes after re-rooting');
$total_length_new = $tree->total_branch_length;
$eps = 0.001 * $total_length_new;    # tolerance for checking length
cmp_ok($total_length_orig, '>=', $tree->total_branch_length - $eps, 'Branch length before rerooting');
cmp_ok($total_length_orig, '<=', $tree->total_branch_length + $eps, 
       'Branch length after rerooting');
is($tree->get_root_node, $a->ancestor->ancestor,'Root is really the ancestor we asked for');

# BFS and DFS search testing
$treeio = Bio::TreeIO->new(-verbose => $verbose,
			     -format => 'newick',
			     -file   => test_input_file('test.nh'));
$tree = $treeio->next_tree;
my ($ct,$n) = (0);
my $let = ord('A');
for $n (  $tree->get_leaf_nodes ) {
    $n->id(chr($let++));
}

for $n ( grep {! $_->is_Leaf } $tree->get_nodes ) {
    $n->id($ct++);
}
# enable for debugging
Bio::TreeIO->new(-format => 'newick')->write_tree($tree) if( $verbose );

my $BFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'b')));
is($BFSorder, '0,1,3,2,C,D,E,F,G,H,A,B', 'BFS traversal order');
my $DFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'd')));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F,G,H', 'DFS travfersal order');


# test some Bio::Tree::TreeFunctionI methods
#find_node tested extensively already
$tree->remove_Node('H');
$DFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'd')));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F,G', 'DFS traversal after removing H');
#get_lineage_nodes tested during get_lca
$tree->splice(-remove_id => 'G');
$DFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'd')));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F', 'DFS traversal after removing G');
$tree->splice(-remove_id => [('E', 'F')], -keep_id => 'F');
$DFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'd')));
is($DFSorder, '0,1,2,A,B,C,D,F', 'DFS traversal after removing F');
$tree->splice(-keep_id => [qw(0 1 2 A B C D)]);
$DFSorder = join(",", map { $_->id } ( $tree->get_nodes(-order => 'd')));
is($DFSorder, '0,1,2,A,B,C,D', 'DFS after removing all but 0,1,2,A,B,C,D');
#get_lca, merge_lineage, contract_linear_paths tested in in Taxonomy.t


# try out the id to bootstrap copy method
$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file   => test_input_file('bootstrap.tre'));
$tree = $treeio->next_tree;
my ($test_node) = $tree->find_node(-id => 'A');
is($test_node->ancestor->id, 90,'Testing bootstrap copy');
is($test_node->ancestor->ancestor->id, '25','Testing bootstrap copy');
$tree->move_id_to_bootstrap;
is($test_node->ancestor->id, '','Testing bootstrap copy');
is($test_node->ancestor->bootstrap, '90', 'Testing bootstrap copy');
is($test_node->ancestor->ancestor->id, '', 'Testing bootstrap copy');
is($test_node->ancestor->ancestor->bootstrap, '25', 'Testing bootstrap copy');

# change TreeIO to parse 
$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file   => test_input_file('bootstrap.tre'),
			   -internal_node_id => 'bootstrap');
$tree = $treeio->next_tree;
($test_node) = $tree->find_node(-id => 'A');
is($test_node->ancestor->id, '','Testing auto-boostrap copy during parse');
is($test_node->ancestor->ancestor->id, '',
   'Testing auto-boostrap copy during parse');
is($test_node->ancestor->bootstrap, '90',
   'Testing auto-boostrap copy during parse');
is($test_node->ancestor->ancestor->bootstrap, '25', 
   'Testing auto-boostrap copy during parse');

__DATA__
(D,(C,(A,B)));
(I,((D,(C,(A,B))),(E,(F,G))));
(((A:0.3,B:2.1):0.45,C:0.7),D:4);
(A:0.031162,((((((B:0.022910,C:0.002796):0.010713,(D:0.015277,E:0.020484):0.005336):0.005588,((F:0.013293,(G:0.018374,H:0.003108):0.005318):0.006047,I:0.014607):0.001677):0.004196,(((((J:0.003307,K:0.001523):0.011884,L:0.006960):0.006514,((M:0.001683,N:0.000100):0.002226,O:0.007085):0.014649):0.008004,P:0.037422):0.005201,(Q:0.000805,R:0.000100):0.015280):0.005736):0.004612,S:0.042283):0.017979,(T:0.006883,U:0.016655):0.040226):0.014239,((((((V:0.000726,W:0.000100):0.028490,((((X:0.011182,Y:0.001407):0.005293,Z:0.011175):0.004701,AA:0.007825):0.016256,BB:0.029618):0.008146):0.004279,CC:0.035012):0.060215,((((((DD:0.014933,(EE:0.008148,FF:0.000100):0.015458):0.003891,GG:0.010996):0.001489,(HH:0.000100,II:0.000100):0.054265):0.003253,JJ:0.019722):0.013796,((KK:0.001960,LL:0.004924):0.013034,MM:0.010071):0.043273):0.011912,(NN:0.031543,OO:0.018307):0.059182):0.026517):0.011087,((PP:0.000100,QQ:0.002916):0.067214,(RR:0.064486,SS:0.013444):0.011613):0.050846):0.015644,((TT:0.000100,UU:0.009287):0.072710,(VV:0.009242,WW:0.009690):0.035346):0.042993):0.060365);
