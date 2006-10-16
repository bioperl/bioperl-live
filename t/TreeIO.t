# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;
use vars qw($NUMTESTS);
use strict;
BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}

	use Test;
	$NUMTESTS = 60;
	plan tests => $NUMTESTS;
}

if( $error == 1 ) {
    exit(0);
}

use vars qw($FILE1 $FILE2 $FILE3);
use File::Spec;
$FILE1= 'testnewick.phylip';
$FILE2= 'testlarge.phy';
$FILE3= File::Spec->catfile(qw(t testsvg.svg));

END {
	unlink $FILE1;
	unlink $FILE2;
	unlink $FILE3;
}
use Bio::TreeIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'newick',
			     -file   => File::Spec->catfile
			     (qw(t data cysprot1b.newick)));

ok($treeio);
my $tree = $treeio->next_tree;
ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

my @nodes = $tree->get_nodes;
ok(@nodes, 6);
my ($rat) = $tree->find_node('CATL_RAT');
ok($rat);
ok($rat->branch_length, '0.12788');
# move the id to the bootstap
 ok($rat->ancestor->bootstrap($rat->ancestor->id), '95');
 $rat->ancestor->id('');
# maybe this can be auto-detected, but then can't distinguish
# between internal node labels and bootstraps...
ok($rat->ancestor->bootstrap, '95');
ok($rat->ancestor->branch_length, '0.18794');
ok($rat->ancestor->id, '');

if($verbose ) {
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
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick',
			  -file   => ">$FILE1");
$treeio->write_tree($tree);
undef $treeio;
ok( -s $FILE1 );
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick',
			  -file   => Bio::Root::IO->catfile('t','data', 
							    'LOAD_Ccd1.dnd'));
ok($treeio);
$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
ok(@nodes, 52);

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

ok($tree->total_branch_length, 7.12148);
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format => 'newick', 
			  -file   => ">$FILE2");
$treeio->write_tree($tree);
undef $treeio;
ok(-s $FILE2);
$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -format  => 'newick',
			  -file    => Bio::Root::IO->catfile('t','data','hs_fugu.newick'));
$tree = $treeio->next_tree();
@nodes = $tree->get_nodes();
ok(@nodes, 5);
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
ok($saw, 1, "Did not see $vals[0] as expected\n");
if( $verbose ) {
	foreach my $node ( @nodes ) {
		print "\t", $node->id, "\n";
	}
}

$treeio = new Bio::TreeIO(-format => 'newick', 
								  -fh => \*DATA);
my $treeout = new Bio::TreeIO(-format => 'tabtree');
my $treeout2 = new Bio::TreeIO(-format => 'newick');

$tree = $treeio->next_tree;

if( $verbose > 0  ) {
    $treeout->write_tree($tree);
    $treeout2->write_tree($tree);
}

$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -file   => Bio::Root::IO->catfile('t','data', 
							    'test.nhx'));

if( eval "require SVG::Graph; 1;" ) {
	my $treeout3 = new Bio::TreeIO(-format => 'svggraph',
											 -file => ">$FILE3");
	ok($treeout3);
	eval {$treeout3->write_tree($tree);};
	ok (-e $FILE3);
} else {
    for ( 1..2 ) {
	skip("skipping SVG::Graph output, SVG::Graph not installed",2);
    }
}

ok($treeio);
$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
ok(@nodes, 13, scalar @nodes);

my $adhy = $tree->find_node('ADHY');
ok($adhy->branch_length, 0.1);
ok(($adhy->get_tag_values('S'))[0], 'nematode');
ok(($adhy->get_tag_values('E'))[0], '1.1.1.1');

# try lintree parsing
$treeio = new Bio::TreeIO(-format => 'lintree',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data crab.njb)));

my (@leaves, $node);
while( $tree = $treeio->next_tree ) {

	ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

	@nodes = $tree->get_nodes;

	@leaves = $tree->get_leaf_nodes;
	ok(@leaves, 13);
	ok(@nodes, 25);
	($node) = $tree->find_node(-id => '18');
	ok($node);
	ok($node->id, '18');
	ok($node->branch_length, '0.030579');
	ok($node->bootstrap, 998);
}

$treeio = new Bio::TreeIO(-format => 'lintree',
			   -file   => Bio::Root::IO->catfile
			   (qw(t data crab.nj)));

$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
@leaves = $tree->get_leaf_nodes;
ok(@leaves, 13);
ok(@nodes, 25);
($node) = $tree->find_node('18');
ok($node->id, '18');
ok($node->branch_length, '0.028117');

($node) = $tree->find_node(-id => 'C-vittat');
ok($node->id, 'C-vittat');
ok($node->branch_length, '0.087619');
ok($node->ancestor->id, '14');

$treeio = new Bio::TreeIO(-format => 'lintree',
			  -file   => Bio::Root::IO->catfile
			  (qw(t data crab.dat.cn)));

$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
@leaves = $tree->get_leaf_nodes;
ok(@leaves, 13, scalar @leaves);

ok(@nodes, 25, scalar @nodes);
($node) = $tree->find_node('18');
ok($node->id, '18');

ok($node->branch_length, '0.029044');

($node) = $tree->find_node(-id => 'C-vittat');
ok($node->id, 'C-vittat');
ok($node->branch_length, '0.097855');
ok($node->ancestor->id, '14');

if( eval "require IO::String; 1;" ) {
# test nexus tree parsing
    $treeio = Bio::TreeIO->new(-format => 'nexus',
			       -file   => Bio::Root::IO->catfile
			       (qw(t data urease.tre.nexus) ));
    
    $tree = $treeio->next_tree;
    ok($tree);
    ok($tree->id, 'PAUP_1');
    ok($tree->get_leaf_nodes, 6);
    ($node) = $tree->find_node(-id => 'Spombe');
    ok($node->branch_length,0.221404);
    
# test nexus MrBayes tree parsing
    $treeio = Bio::TreeIO->new(-format => 'nexus',
			       -file   => Bio::Root::IO->catfile
			       (qw(t data adh.mb_tree.nexus) ));
    
    $tree = $treeio->next_tree;
    ok($tree);
    ok($tree->id, 'rep.1');
    ok($tree->get_leaf_nodes, 54);
    ($node) = $tree->find_node(-id => 'd.madeirensis');
    ok($node->branch_length,0.039223);
} else{
    for ( 1..8 ) {
	skip("skipping nexus tree parsing, IO::String not installed",1);
    }
}

# bug #1854
# process no-newlined tree
$treeio = Bio::TreeIO->new(-format => 'nexus',
			   -file   => Bio::Root::IO->catfile
			   (qw(t data tree_nonewline.nexus) ));

$tree = $treeio->next_tree;
ok($tree);
ok($tree->find_node('TRXHomo'));


# parse trees with scores

$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file   => Bio::Root::IO->catfile
			   (qw(t data puzzle.tre)));
$tree = $treeio->next_tree;
ok($tree);
ok($tree->score, '-2673.059726');
							     
__DATA__
(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);
