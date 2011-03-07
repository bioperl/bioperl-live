use Bio::Tree::Tree;
my $tree;
 
# Step 1: Load a tree
#   from a file
#$tree = Bio::Tree::Tree->from_file("my-tree.xml"); 
#   from a string
$tree = Bio::Tree::Tree->from_string("(a,(b,c));"); 
#   using the TreeIO system (the format can be autodetected if the -format argument is missing)
#my $treeio = Bio::TreeIO->new(-file => "my-tree.xml", -format => 'phyloxml'); 
#$tree = $treeio->next_tree;
 
# Step 2: Save a tree
#   you can get the Newick string directly
print $tree->newick . "\n";
#   or use a TreeIO writer
$treeio = Bio::TreeIO->new(-file => ">tree-out.xml", -format => 'phyloxml');
$treeio->write_tree($tree);
 
# Step 3: Analyze a tree
print "Leaf count: " . scalar($tree->leaves) . "\n";
print "Node count: " . scalar($tree->nodes) . "\n";
print "Total branch length: " . $tree->total_branch_length . "\n";
print "Max root-to-tip branch length: " . $tree->max_distance_to_leaf . "\n";
print "Max root-to-tip node depth: " . $tree->max_depth_to_leaf . "\n";
 
# Step 4: Modify a tree
#   print a human-readable ASCII diagram
print "Original tree:\n" . $tree->ascii; 
my $orig_root = $tree->root;
#   re-root the tree on the node labeled 'b'
$tree->reroot($tree->find('b'));
print "Rerooted on b:\n" . $tree->ascii;
#   re-root halfway along the branch leading to 'c'
#   (this can be more intuitive, but it adds a new internal
#   node to the tree)
$tree->reroot_above($tree->find('c'), 0.5); 
print "Rerooted above c:\n" . $tree->ascii;
#   return to the old root and remove the internal node created
#   in the previous re-rooting.
$tree->reroot($orig_root);
$tree->contract_linear_paths;

#   use key-value mappings to translate the tree's node labels
my $id_map = {
  'a' => 'Aardvark',
  'b' => 'Banana',
  'c' => 'Coyote'
};
$tree->translate_ids($id_map);
print "Translated IDs:\n" . $tree->ascii;
 
#   add a new node to the tree
#   ... first create a new object of the same class as the root
my $root_node = $tree->root;
my $new_node = new $root_node;
#   ... then add it as third child of the node parental to Banana
$tree->find('Banana')->parent->add_child($new_node);
$new_node->branch_length(1);
$new_node->id('z');
print "New node added:\n" . $tree->ascii;

#   now the tree has a multifurcation -- ask BioPerl to randomly
#   resolve the multifurcation
$tree->force_binary;
print "Forced to binary structure:\n" . $tree->ascii;

$tree = Bio::Tree::Tree->from_string("(a,(b,(c,(d,(e,(f,(g,(h,i))))))));");
print "Full tree:\n" . $tree->ascii;
my $slice = $tree->slice($tree->find('a'), $tree->find('d'), $tree->find('i'));
print "Slice:\n" . $slice->ascii;
my $slice = $tree->slice_by_ids('a', 'c', 'e', 'i');
print "Slice:\n" . $slice->ascii;

