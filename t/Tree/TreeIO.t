# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 78);
    
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

ok my $treeio = Bio::TreeIO->new(-verbose => $verbose,
                 -format => 'newick',
                 -file   => test_input_file('cysprot1b.newick'));

my $tree = $treeio->next_tree;
isa_ok($tree, 'Bio::Tree::TreeI');

my @nodes = $tree->get_nodes;
is(@nodes, 6);
my ($rat) = $tree->find_node('CATL_RAT');
ok($rat);
is($rat->branch_length, '0.12788');
# move the id to the bootstap
is($rat->ancestor->bootstrap($rat->ancestor->id), '95');
$rat->ancestor->id('');
# maybe this can be auto-detected, but then can't distinguish
# between internal node labels and bootstraps...
is($rat->ancestor->bootstrap, '95');
is($rat->ancestor->branch_length, '0.18794');
is($rat->ancestor->id, '');

if ($verbose) {
    foreach my $node ( $tree->get_root_node()->each_Descendent() ) {
        print "node: ", $node->to_string(), "\n";
        my @ch = $node->each_Descendent();
        if( @ch ) {
            print "\tchildren are: \n";
            foreach my $node ( $node->each_Descendent() ) {
                print "\t\t ", $node->to_string(), "\n";
            }
        }
    }
}

my $FILE1 = test_output_file();
$treeio = Bio::TreeIO->new(-verbose => $verbose,
              -format => 'newick',
              -file   => ">$FILE1");
$treeio->write_tree($tree);
undef $treeio;
ok( -s $FILE1 );
$treeio = Bio::TreeIO->new(-verbose => $verbose,
              -format => 'newick',
              -file   => test_input_file('LOAD_Ccd1.dnd'));
ok($treeio);
$tree = $treeio->next_tree;

isa_ok($tree,'Bio::Tree::TreeI');

@nodes = $tree->get_nodes;
is(@nodes, 52);

if( $verbose ) { 
    foreach my $node ( @nodes ) {
        print "node: ", $node->to_string(), "\n";
        my @ch = $node->each_Descendent();
        if( @ch ) {
            print "\tchildren are: \n";
            foreach my $node ( $node->each_Descendent() ) {
                print "\t\t ", $node->to_string(), "\n";
            }
        }
    }
}

is($tree->total_branch_length, 7.12148);
my $FILE2 = test_output_file();
$treeio = Bio::TreeIO->new(-verbose => $verbose,
              -format => 'newick', 
              -file   => ">$FILE2");
$treeio->write_tree($tree);
undef $treeio;
ok(-s $FILE2);
$treeio = Bio::TreeIO->new(-verbose => $verbose,
              -format  => 'newick',
              -file    => test_input_file('hs_fugu.newick'));
$tree = $treeio->next_tree();
@nodes = $tree->get_nodes();
is(@nodes, 5);
# no relable order for the bottom nodes because they have no branchlen
my @vals = qw(SINFRUP0000006110);
my $saw = 0;
foreach my $node ( $tree->get_root_node()->each_Descendent() ) {
    foreach my $v ( @vals ) {
       if( defined $node->id && 
           $node->id eq $v ){ $saw = 1; last; }
    }
    last if $saw;
}
is($saw, 1, "Saw $vals[0] as expected");
if( $verbose ) {
    foreach my $node ( @nodes ) {
        print "\t", $node->id, "\n" if $node->id;
    }
}

$treeio = Bio::TreeIO->new(-format => 'newick', 
                                  -fh => \*DATA);
my $treeout = Bio::TreeIO->new(-format => 'tabtree');
my $treeout2 = Bio::TreeIO->new(-format => 'newick');

$tree = $treeio->next_tree;

if( $verbose > 0  ) {
    $treeout->write_tree($tree);
    $treeout2->write_tree($tree);
}

$treeio = Bio::TreeIO->new(-verbose => $verbose,
              -file   => test_input_file('test.nhx'));

SKIP: {
    test_skip(-tests => 2, -requires_module => 'SVG::Graph');
    my $FILE3 = test_output_file();
    my $treeout3 = Bio::TreeIO->new(-format => 'svggraph',
                                             -file => ">$FILE3");
    ok($treeout3);
    eval {$treeout3->write_tree($tree);};
    ok (-s $FILE3);
}

ok($treeio);
$tree = $treeio->next_tree;

isa_ok($tree, 'Bio::Tree::TreeI');

@nodes = $tree->get_nodes;
is(@nodes, 12, "Total Nodes");

my $adhy = $tree->find_node('ADHY');
is($adhy->branch_length, 0.1);
is(($adhy->get_tag_values('S'))[0], 'nematode');
is(($adhy->get_tag_values('E'))[0], '1.1.1.1');

# try lintree parsing
$treeio = Bio::TreeIO->new(-format => 'lintree',
                  -file   => test_input_file('crab.njb'));

my (@leaves, $node);
while( $tree = $treeio->next_tree ) {

    isa_ok($tree, 'Bio::Tree::TreeI');

    @nodes = $tree->get_nodes;

    @leaves = $tree->get_leaf_nodes;
    is(@leaves, 13);
#/maj   is(@nodes, 25);
    is(@nodes, 24); # this is clear from the datafile and counting \maj
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

@nodes = $tree->get_nodes;
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
is(@nodes, 24, "All nodes"); 
($node) = $tree->find_node('18');
is($node->id, '18');

is($node->branch_length, '0.029044');

($node) = $tree->find_node(-id => 'C-vittat');
is($node->id, 'C-vittat');
is($node->branch_length, '0.097855');
is($node->ancestor->id, '14');

SKIP: {
    test_skip(-tests => 8, -requires_module => 'IO::String');
    
    # test nexus tree parsing
    $treeio = Bio::TreeIO->new(-format => 'nexus',
                               -verbose => $verbose,
                   -file   => test_input_file('urease.tre.nexus'));
    
    $tree = $treeio->next_tree;
    ok($tree);
    is($tree->id, 'PAUP_1');
    is($tree->get_leaf_nodes, 6);
    ($node) = $tree->find_node(-id => 'Spombe');
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
$treeio = Bio::TreeIO->new(-format => 'nexus',
                           -verbose => $verbose,
               -file   => test_input_file('tree_nonewline.nexus'));

$tree = $treeio->next_tree;
ok($tree);
ok($tree->find_node('TRXHomo'));


# parse trees with scores

$treeio = Bio::TreeIO->new(-format => 'newick',
               -file   => test_input_file('puzzle.tre'));
$tree = $treeio->next_tree;
ok($tree);
is($tree->score, '-2673.059726');

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

# bug #2869
#

# proper way (Tree isn't GC'd)
$tree = Bio::TreeIO->new(-format => 'newick',
               -verbose => $verbose,
               -file   => test_input_file('bug2869.tree'))->next_tree;

my $root = $tree->get_root_node;

isa_ok($root, 'Bio::Tree::NodeI');

my $total1 = 0;
for my $child ($root->get_Descendents) {
    $total1++;
}

is($total1, 198);

undef $tree; # GC

$root = Bio::TreeIO->new(-format => 'newick',
               -verbose => $verbose,
               -file   => test_input_file('bug2869.tree'),
               -no_cleanup => 1)->next_tree->get_root_node;

isa_ok($root, 'Bio::Tree::NodeI');

TODO: {
    local $TODO = 'The nodes are garbage-collected away b/c Tree isn\'t retained in memory';
    my $total2 = 0;
    for my $child ($root->get_Descendents) {
        $total2++;
    }
    
    is($total2, $total1);
}

__DATA__
(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);
