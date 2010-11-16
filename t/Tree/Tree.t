# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin();
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

my $treeio = Bio::TreeIO->new(-verbose => $verbose,
			     -format => 'nhx',
			     -file   => test_input_file('test.nhx'));
my $tree = $treeio->next_tree;

# tests for tags
ok ! $tree->has_tag('test');
is $tree->add_tag_value('test','a'), 1;
ok $tree->has_tag('test');
is $tree->add_tag_value('test','b'), 2;
my @tags = $tree->get_tag_values('test');
is scalar @tags, 2;
is scalar $tree->get_tag_values('test'), 'a', 'retrieve the first value';
is $tree->remove_tag('test2'), 0;
is $tree->remove_tag('test'), 1;
ok ! $tree->has_tag('test');
is $tree->set_tag_value('test',('a','b','c')), 3;
is $tree->remove_all_tags(), undef;
ok ! $tree->has_tag('test');


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

is($tree->is_monophyletic(\@hADH,$n4),1,'Test Monophyly');

my @mixgroup = ( $tree->find_node('hADH1'),
		 $tree->find_node('yADH2'),
		 $tree->find_node('yADH3'),
		 );

my ($iADHX) = $tree->find_node('iADHX');

# test height
is($iADHX->height, 0,'Height');
is($iADHX->distance_to_root,0.22,'Distance to root');
isnt( $tree->is_monophyletic(\@mixgroup,$iADHX),1, 'non-monophyletic group');

# binary tree?
is $tree->is_binary, 0, 'not a binary tree';
#print STDERR $tree->ascii;
is scalar $tree->nodes, 12, '12 nodes';
$tree->verbose(-1);
#print STDERR $tree->ascii;
$tree->force_binary;
#print STDERR $tree->ascii;
is $tree->is_binary, 1, 'after force_binary() it is';
isnt(scalar $tree->nodes, 12, 'and there are more nodes (not 12 anymore)');

my $in = Bio::TreeIO->new(-format => 'newick',
			 -fh     => \*DATA);
$tree = $in->next_tree;
my ($a,$b,$c,$d) = ( $tree->find_node('A'),
		     $tree->find_node('B'),
		     $tree->find_node('C'),
		     $tree->find_node('D'));

is($tree->is_monophyletic([$b,$c],$d),0, 'B,C are NOT Monophyletic w/ D');

is($tree->is_monophyletic([$b,$a],$c),1,'A,B are Monophyletic with C');

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
is( $tree->is_monophyletic([$b,$f],$d),0,'B,F are not Monophyletic' );
is($tree->is_monophyletic([$b,$a],$f),1, 'A,B are Monophyletic');

# tests for paraphyly
is(  $tree->is_paraphyletic([$a,$b,$c],$i), 1, 'A,B,C are paraphyletic w/ I');
is(  $tree->is_paraphyletic([$a,$b],$d), 1, 'A,B are paraphyletic w/ D');
is(  $tree->is_paraphyletic([$a,$b],$i), 1, 'A,B are paraphyletic w/ I');
is(  $tree->is_paraphyletic([$a,$b],$c), 1, 'A,B are paraphyletic w/ C');
is(  $tree->is_paraphyletic([$a,$f,$e],$i), 1, 'A,F,E are paraphyletic w/ I');
is(  $tree->is_paraphyletic([$a,$f,$e],$d), 0,'A,F,E are polyphyletic w/ D, so not paraphyletic');

# tests for polyphyly
is(  $tree->is_polyphyletic([$a,$f,$e],$d), 1,'A,F,E are polyphyletic w/ D');

# Test cloning.
my $tree_two = $tree->clone;
is($tree->as_text('newick'),$tree_two->as_text('newick'),"Cloned tree equals orig");

# test for rerooting the tree
my $out = Bio::TreeIO->new(-format => 'newick', 
			   -fh => \*STDERR, 
			   -noclose => 1);
$tree = $in->next_tree;
$tree->verbose( -1 ) unless $verbose;
my $node_cnt_orig = scalar($tree->nodes);
# reroot on an internal node: should work fine
$a = $tree->find_node('A');
# removing node_count checks because re-rooting can change the
# number of internal nodes (if it is done correctly)
my $total_length_orig = $tree->total_branch_length;
is $tree->total_branch_length, $tree->subtree_length, 
    "subtree_length() without attributes is an alias to total_branch_lenght()";
cmp_ok($total_length_orig, '>',$a->ancestor->subtree_length, 
       'Length of the tree is larger than length of a subtree');
$out->write_tree($tree) if $verbose;
is($tree->reroot($a),1, 'Can re-root with A as outgroup');
$out->write_tree($tree) if $verbose;
is($node_cnt_orig, scalar($tree->nodes), 'Count the number of nodes');
my $total_length_new = $tree->total_branch_length;
my $eps = 0.001 * $total_length_new;	# tolerance for checking length
warn("orig total len ", $total_length_orig, "\n") if $verbose;
warn("new  total len ", $tree->total_branch_length,"\n") if $verbose;
# according to retree in phylip these branch lengths actually get larger
# go figure...
# this should be fixed now/maj
ok(($total_length_orig >= $tree->total_branch_length - $eps) &&
   ($total_length_orig <= $tree->total_branch_length + $eps),'same length');

# prob with below: rerooted tree on node A at line 146; so $a IS root
#/maj is($tree->get_root_node, $a->ancestor, "Root node is A's ancestor");
is($tree->get_root_node, $a, "Root node is A");

# former test expected the old behavior of reroot; here is the new
# test/maj
my $desc = ($a->each_Descendent)[0];
$tree->reroot_above($desc,0.5);
#$tree->reroot($newroot);
is($tree->get_root_node, $a->ancestor, "Root node is A's ancestor");

# try to reroot on an internal, will result in there being 1 less node
# Rerooting should be an invariant operation with respect to node number!/maj
# the test show that it now is, because the secret removal of nodes 
# no longer occurs

$a = $tree->find_node('C')->ancestor;
$out->write_tree($tree) if $verbose;
is($tree->reroot($a),1, "Can reroot with C's ancsestor");
$out->write_tree($tree) if $verbose;
#/maj is($node_cnt_orig, scalar($tree->get_nodes), 'Check to see that node count is correct after an internal node was removed after this re-rooting');
# but we did add a new node at line 166, so
is($node_cnt_orig+1, scalar($tree->get_nodes), 'Node count correct');
warn("orig total len ", $total_length_orig, "\n") if $verbose;
warn("new  total len ", $tree->total_branch_length,"\n") if $verbose;
cmp_ok($total_length_orig, '>=', $tree->total_branch_length - $eps, 
       'Total original branch length is what it is supposed to be');
# branch length should also be invariant w/r to rerooting...
cmp_ok($total_length_orig, '<=',$tree->total_branch_length + $eps, 
       'Updated total branch length after the reroot');
# again, we rerooted ON THE NODE, so $a IS the root./maj
is($tree->get_root_node, $a, 'Make sure root is really what we asked for');

# try to reroot on new root: should fail
#/maj  $a = $tree->get_root_node;
isnt( $tree->reroot($a),1, 'Testing for failed re-rerooting');

# try a more realistic tree
$tree = $in->next_tree;
$a = $tree->find_node('VV');
$node_cnt_orig = scalar($tree->get_nodes);
$total_length_orig = $tree->total_branch_length;
$out->write_tree($tree) if $verbose;
is($tree->reroot($a),1, 'Test that rooting succeeded'); #mod /maj
$out->write_tree($tree) if $verbose;
# node number should be invariant after reroot/maj
is($node_cnt_orig, scalar($tree->get_nodes), 'Test that re-rooted tree has proper number of nodes after re-rooting'); #mod /maj
$total_length_new = $tree->total_branch_length;
$eps = 0.001 * $total_length_new;    # tolerance for checking length
cmp_ok($total_length_orig, '>=', $tree->total_branch_length - $eps, 'Branch length before rerooting');
cmp_ok($total_length_orig, '<=', $tree->total_branch_length + $eps, 
       'Branch length after rerooting');
is($tree->get_root_node, $a,'Root is really the ancestor we asked for'); #mod /maj

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


my $BFSorder = join(",", map { $_->id } ( $tree->nodes_breadth_first));
is($BFSorder, '0,1,3,2,C,D,E,F,G,H,A,B', 'BFS traversal order');
my $DFSorder = join(",", map { $_->id } ( $tree->nodes_depth_first));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F,G,H', 'DFS travfersal order');


# test some Bio::Tree::TreeFunctionI methods
#find_node tested extensively already
$tree->remove_Node('H');
$DFSorder = join(",", map { $_->id } ( $tree->nodes_depth_first));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F,G', 'DFS traversal after removing H');
#get_lineage_nodes tested during get_lca
#$tree->splice(-remove_id => 'G');
$tree->find('G')->splice();
$DFSorder = join(",", map { $_->id } ( $tree->nodes_depth_first));
is($DFSorder, '0,1,2,A,B,C,D,3,E,F', 'DFS traversal after removing G');
$tree->find('E')->splice();
#$tree->splice(-remove_id => [('E', 'F')], -keep_id => 'F');
$DFSorder = join(",", map { $_->id } ( $tree->nodes_depth_first));
# the node '3' is not explicitly removed, so it should still be there
# I suspect that it disappeared before was due to the previously
# automatic removal of internal degree 2 nodes../maj
is($DFSorder, '0,1,2,A,B,C,D,3,F', 'DFS traversal after removing E');
$tree->find('3')->splice();
$tree->find('F')->splice();
#$tree->splice(-keep_id => [qw(0 1 2 A B C D)]);
$DFSorder = join(",", map { $_->id } ( $tree->nodes_depth_first));
is($DFSorder, '0,1,2,A,B,C,D', 'DFS after removing all but 0,1,2,A,B,C,D');

my $lca = $tree->find('A')->lca($tree->find('D'));
is($lca->id,'1',"LCA of A and D");

#get_lca, merge_lineage, contract_linear_paths tested in in Taxonomy.t


# try out the id to bootstrap copy method
$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file   => test_input_file('bootstrap.tre'));
$tree = $treeio->next_tree;
my ($test_node) = $tree->root->find_by_id('A');
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
($test_node) = $tree->root->find_by_id('A');
is($test_node->ancestor->id, '','Testing auto-boostrap copy during parse');
is($test_node->ancestor->ancestor->id, '',
   'Testing auto-boostrap copy during parse');
is($test_node->ancestor->bootstrap, '90',
   'Testing auto-boostrap copy during parse');
is($test_node->ancestor->ancestor->bootstrap, '25', 
   'Testing auto-boostrap copy during parse');

sub tree_from_string {
    my $string = shift;
    return Bio::TreeIO->new(-string=>$string)->next_tree;
}

sub test_flip {
    my $string = shift;
    my $flipped_string = shift;
    my $tree = tree_from_string($string);
    $tree->get_root_node->flip_subtree;
    my $test = $tree->as_text('newick');
    is($test,$flipped_string,"Testing flip [$string] should flip to [$flipped_string]");
}
test_flip("(a,(b,c));","((c,b),a);");
test_flip("(a,b,c,d,e,f,g);","(g,f,e,d,c,b,a);");

sub test_rev {
    my $string = shift;
    my $flipped_string = shift;
    my $tree = tree_from_string($string);
    $tree->get_root_node->reverse_children;
    my $test = $tree->as_text('newick');
    is($test,$flipped_string,"Testing reverse [$string] should reverse to [$flipped_string]");
}
test_rev("(a,b);","(b,a);");
test_rev("(a,(b,c,d));","((b,c,d),a);");


# Do some re-rooting testing.
sub test_reroot {
    my $tree = shift;
    my $node = shift;

    my $orig_str = $tree->as_text('newick');
    my $orig_root = $tree->root;
    $tree->reroot($node);
    $tree->reroot($orig_root);
    my $str = $tree->as_text('newick');
    is($orig_str,$str,"Testing rerooting tree on node ".$node->id);
}
$tree = tree_from_string("(a,(b,(c,(d,(e,(f,(g,h)))))));");
# Test a re-root roundtrip for each of the tree's nodes.
foreach my $node ($tree->root->nodes) {
  test_reroot($tree,$node);
}

$tree = tree_from_string("(a[&&NHX:type=whatever],(b[&&NHX:type=however],(c[&&NHX:type=evermore],(d,(e,(f,(g,h)))))));");
my @found;
@found = $tree->root->find_by_tag_value('type','whatever');
is(scalar(@found),1,"Find node by specific tag value");
@found = $tree->root->find_by_tag_value('type','^.*ever$');
is(scalar(@found),2,"Find nodes by tag_value regex");
@found = $tree->root->find_by_tag_value('type','ever');
is(scalar(@found),3,"Find nodes by tag_value regex");


# Test force binary script.
$tree = tree_from_string("(a,b,c,d,e,f,g,h);");
#print STDERR $tree->ascii;
$tree->root->force_binary();
#print STDERR $tree->ascii;

$tree = tree_from_string("(((((((a,b)c)d)e)f)g)h)i;");

$tree->root->remove_elbow_nodes;
is($tree->as_text,'((a,b)c)i;',"Keep root node after removing linear paths");

$lca = $tree->find('a')->lca($tree->find('b'));
is($lca->id,'c',"LCA after contracting linear paths");

$tree = tree_from_string("((a,b),(c,d))e;");
$lca = $tree->get_lca($tree->nodes);
is($lca,$tree->root,"LCA of all nodes is the root node");


# Distance measurement tests
$tree = tree_from_string("(a:1,b:2,c:3,d:4,e:5)root:2;");
is($tree->distance($tree->find('a'),$tree->find('b')),3,"Distance betwen a and b");
is($tree->distance($tree->find('a'),$tree->find('e')),6,"Distance betwen a and e");
is($tree->distance($tree->find('d'),$tree->find('e')),9,"Distance betwen d and e");



done_testing();

__DATA__
(D,(C,(A,B)));
(I,((D,(C,(A,B)x)y),(E,(F,G))));
(((A:0.3,B:2.1):0.45,C:0.7),D:4);
(A:0.031162,((((((B:0.022910,C:0.002796):0.010713,(D:0.015277,E:0.020484):0.005336):0.005588,((F:0.013293,(G:0.018374,H:0.003108):0.005318):0.006047,I:0.014607):0.001677):0.004196,(((((J:0.003307,K:0.001523):0.011884,L:0.006960):0.006514,((M:0.001683,N:0.000100):0.002226,O:0.007085):0.014649):0.008004,P:0.037422):0.005201,(Q:0.000805,R:0.000100):0.015280):0.005736):0.004612,S:0.042283):0.017979,(T:0.006883,U:0.016655):0.040226):0.014239,((((((V:0.000726,W:0.000100):0.028490,((((X:0.011182,Y:0.001407):0.005293,Z:0.011175):0.004701,AA:0.007825):0.016256,BB:0.029618):0.008146):0.004279,CC:0.035012):0.060215,((((((DD:0.014933,(EE:0.008148,FF:0.000100):0.015458):0.003891,GG:0.010996):0.001489,(HH:0.000100,II:0.000100):0.054265):0.003253,JJ:0.019722):0.013796,((KK:0.001960,LL:0.004924):0.013034,MM:0.010071):0.043273):0.011912,(NN:0.031543,OO:0.018307):0.059182):0.026517):0.011087,((PP:0.000100,QQ:0.002916):0.067214,(RR:0.064486,SS:0.013444):0.011613):0.050846):0.015644,((TT:0.000100,UU:0.009287):0.072710,(VV:0.009242,WW:0.009690):0.035346):0.042993):0.060365);
