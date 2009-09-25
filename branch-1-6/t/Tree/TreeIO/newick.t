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

# parse trees with scores

$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file   => test_input_file('puzzle.tre'));
$tree = $treeio->next_tree;
ok($tree);
is($tree->score, '-2673.059726');


# no semi-colon

$treeio = Bio::TreeIO->new(-format => 'newick', 
			   -file=> test_input_file('semicolon.newick'));
$tree = $treeio->next_tree;
ok($tree);
is($tree->get_nodes, 15);

$treeio = Bio::TreeIO->new(-format => 'newick', 
			   -file=> test_input_file('no_semicolon.newick'));
$tree = $treeio->next_tree;
ok($tree);
is($tree->get_nodes, 15);
