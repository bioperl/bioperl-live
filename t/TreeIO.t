# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
    plan tests => 17; 

#    eval { require XML::Parser::PerlSAX; };
#    if( $@ ) {
#	print STDERR "XML::Parser::PerlSAX not loaded. This means TreeIO::phyloxml test cannot be executed. Skipping\n";
#	foreach ( 1..43 ) {
#	    skip(1,1);
#	}
#       $error = 1;
#	
#    } 

}

if( $error == 1 ) {
    exit(0);
}

use vars qw($FILE1 $FILE2);

$FILE1= 'testnewick.phylip';
$FILE2= 'testlarge.phy';

END { 
	unlink $FILE1;
	unlink $FILE2;
}
use Bio::TreeIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'} || 0;
ok(1);

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'newick',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'cysprot1b.newick'));

ok($treeio);
my $tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

my @nodes = $tree->get_nodes;
ok(@nodes, 6);

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

if($verbose ) { 
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
	   if( $node->id eq $v ){ $saw = 1; last; }
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
@nodes = $tree->get_nodes;

my( $i, $c, $g);

for ($i = 0; $i <= $#nodes; $i++) {
    next unless defined $nodes[$i]->id;
    if ($nodes[$i]->id eq 'C') {
	$c = $i;
    }
    if ($nodes[$i]->id eq 'G') {
	$g = $i;
    }
}
$nodes[$c]->ancestor;
$nodes[$g]->ancestor;
my $cancestor = $nodes[$c]->ancestor;
my $gancestor = $nodes[$g]->ancestor; 
$cancestor->id('C-ancestor'); # let's provide a way to test if we suceeded
$gancestor->id('G-ancestor'); # in our swapping

$cancestor->remove_Descendent($nodes[$c]);
$gancestor->remove_Descendent($nodes[$g]);
$cancestor->add_Descendent($nodes[$g],1);
$gancestor->add_Descendent($nodes[$c],1);

@nodes = $tree->get_nodes();

for ($i = 0; $i <= $#nodes; $i++) {
    next unless defined $nodes[$i]->id;
    if ($nodes[$i]->id eq 'C') {
	ok($nodes[$i]->ancestor->id, 'G-ancestor');
	$c = $i;
    }
    if ($nodes[$i]->id eq 'G') {
	$g = $i;
	ok($nodes[$i]->ancestor->id, 'C-ancestor');
    }
}

if( $verbose > 0  ) {
    $treeout2->write_tree($tree);
}

$treeio = new Bio::TreeIO(-verbose => $verbose,
			  -file   => Bio::Root::IO->catfile('t','data', 
							    'test.nhx'));

ok($treeio);
$tree = $treeio->next_tree;

ok(ref($tree) && $tree->isa('Bio::Tree::TreeI'));

@nodes = $tree->get_nodes;
ok(@nodes, 12, scalar @nodes);



__DATA__
(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);
