# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    $error = 0; 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 22;
    plan tests => $TESTCOUNT;
}

use Bio::TreeIO;
my $verbose = 0;

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'nhx',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'test.nhx'));
my $tree = $treeio->next_tree;

my @nodes = $tree->find_node('ADH2');
ok(@nodes, 2);

if( $verbose ) {
    $treeio = new Bio::TreeIO(-verbose => $verbose,
			      -format => 'nhx',
			      );
    $treeio->write_tree($tree);
    print "nodes are: \n",
    join(", ", map {  $_->id . ":". (defined $_->branch_length ? 
				     $_->branch_length : '' ) } @nodes), "\n";
}

$treeio = new Bio::TreeIO(-format => 'newick',
			  -file   => Bio::Root::IO->catfile('t','data',
							    'test.nh'));
$tree = $treeio->next_tree;
if( $verbose ) { 
    my $out = new Bio::TreeIO(-format => 'tabtree');
    
    $out->write_tree($tree);
}

my @hADH = ( $tree->find_node('hADH1'),
	     $tree->find_node('hADH2') );
my ($n4) = $tree->find_node('yADH4');

ok($tree->is_monophyletic(-nodes    => \@hADH,
			  -outgroup => $n4));

my @mixgroup = ( $tree->find_node('hADH1'),
		 $tree->find_node('yADH2'),
		 $tree->find_node('yADH3'),
		 );

my ($iADHX) = $tree->find_node('iADHX');

ok(! $tree->is_monophyletic(-nodes   => \@mixgroup,
			    -outgroup=> $iADHX));

my $in = new Bio::TreeIO(-format => 'newick',
			 -fh     => \*DATA);
$tree = $in->next_tree;
my ($a,$b,$c,$d) = ( $tree->find_node('A'),
			   $tree->find_node('B'),
			   $tree->find_node('C'),
			   $tree->find_node('D'));

ok($tree->is_monophyletic(-nodes => [$b,$c],
			  -outgroup => $d));

ok($tree->is_monophyletic(-nodes => [$b,$a],
			  -outgroup => $d) );

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
ok(! $tree->is_monophyletic(-nodes => [$b,$f],
			    -outgroup => $d) );

ok($tree->is_monophyletic(-nodes => [$b,$a],
			  -outgroup => $f));

# test for paraphyly

ok(  $tree->is_paraphyletic(-nodes => [$a,$b,$c],
			   -outgroup => $d), 0);

ok(  $tree->is_paraphyletic(-nodes => [$a,$f,$e],
			   -outgroup => $i), 1);


# test for rerooting the tree

$tree = $in->next_tree;
$tree->verbose(-1) unless $DEBUG;
# reroot on an internal node: should work fine
$a = $tree->find_node('A');
my $node_cnt_orig = scalar($tree->get_nodes);
my $total_length_orig = $tree->total_branch_length;
ok($tree->reroot($a->ancestor) eq '1');
ok($node_cnt_orig, scalar($tree->get_nodes));
my $total_length_new = $tree->total_branch_length;
my $eps = 0.001 * $total_length_new;	# tolerance for checking length
ok(($total_length_orig >= $tree->total_branch_length - $eps)
   and ($total_length_orig <= $tree->total_branch_length + $eps));
ok($tree->get_root_node, $a->ancestor);
# try to reroot on a leaf: should end up rerooting on its ancestor
$a = $tree->find_node('C');
ok($tree->reroot($a) eq '1');
ok($node_cnt_orig, scalar($tree->get_nodes));
ok(($total_length_orig >= $tree->total_branch_length - $eps)
   and ($total_length_orig <= $tree->total_branch_length + $eps));
ok($tree->get_root_node, $a->ancestor);
# try to reroot on existing root: should fail
$a = $tree->get_root_node;
ok($tree->reroot($a) eq '0');
# try a more realistic tree
$tree = $in->next_tree;
$a = $tree->find_node('VV');
$node_cnt_orig = scalar($tree->get_nodes);
$total_length_orig = $tree->total_branch_length;
ok($tree->reroot($a->ancestor) eq '1');
ok($node_cnt_orig, scalar($tree->get_nodes));
$total_length_new = $tree->total_branch_length;
$eps = 0.001 * $total_length_new;    # tolerance for checking length
ok(($total_length_orig >= $tree->total_branch_length - $eps)
   and ($total_length_orig <= $tree->total_branch_length + $eps));
ok($tree->get_root_node, $a->ancestor);

__DATA__
(D,(C,(A,B)));
(I,((D,(C,(A,B))),(E,(F,G))));
((A:0.3,B:2.1):0.000003,C:1e-2,D:4);
(A:0.031162,((((((B:0.022910,C:0.002796):0.010713,(D:0.015277,E:0.020484):0.005336):0.005588,((F:0.013293,(G:0.018374,H:0.003108):0.005318):0.006047,I:0.014607):0.001677):0.004196,(((((J:0.003307,K:0.001523):0.011884,L:0.006960):0.006514,((M:0.001683,N:0.000100):0.002226,O:0.007085):0.014649):0.008004,P:0.037422):0.005201,(Q:0.000805,R:0.000100):0.015280):0.005736):0.004612,S:0.042283):0.017979,(T:0.006883,U:0.016655):0.040226):0.014239,((((((V:0.000726,W:0.000100):0.028490,((((X:0.011182,Y:0.001407):0.005293,Z:0.011175):0.004701,AA:0.007825):0.016256,BB:0.029618):0.008146):0.004279,CC:0.035012):0.060215,((((((DD:0.014933,(EE:0.008148,FF:0.000100):0.015458):0.003891,GG:0.010996):0.001489,(HH:0.000100,II:0.000100):0.054265):0.003253,JJ:0.019722):0.013796,((KK:0.001960,LL:0.004924):0.013034,MM:0.010071):0.043273):0.011912,(NN:0.031543,OO:0.018307):0.059182):0.026517):0.011087,((PP:0.000100,QQ:0.002916):0.067214,(RR:0.064486,SS:0.013444):0.011613):0.050846):0.015644,((TT:0.000100,UU:0.009287):0.072710,(VV:0.009242,WW:0.009690):0.035346):0.042993):0.060365);
