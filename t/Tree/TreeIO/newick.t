# -*-Perl-*- Test Harness script for Bioperl
# $Id: TreeIO.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
    
  test_begin(-tests => 51);

  use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

ok my $treeio = Bio::TreeIO->new(-verbose => $verbose,
                                 -format => 'newick',
                                 -file => test_input_file('cysprot1b.newick'));

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

# initial AMPHORA2 tests
$treeio = Bio::TreeIO->new(-format => 'newick',
			   -file=> test_input_file('amphora.newick'));
$tree = $treeio->next_tree;
ok($tree);
is($tree->get_nodes, 5);

test_roundtrip('((a,b),c);','Round trip: simple newick');
test_roundtrip('(a:1,b:2,c:3,d:4)TEST:1.2345;','Round trip: Root node branch length');
test_roundtrip('(a:1,b:2,c:3,d:4):1.2345;','Round trip: Root node branch length');
test_roundtrip('(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;','Round trip: from Wikipedia');
test_roundtrip('(a:1,b:2):0.0;','Branch length on root');
test_roundtrip('(a:1,b:2):0.001;','Tiny branch length on root');
test_roundtrip('(a:0,b:00):0.0;','Zero branch lenghts');

# From Wikipedia:
test_roundtrip('(,,(,));','wkp blank tree');
test_roundtrip('(A,B,(C,D));','wkp only leaves labeled');
test_roundtrip('(A,B,(C,D)E)F;','wkp all nodes labeled');
test_roundtrip('(:0.1,:0.2,(:0.3,:0.4):0.5);','wkp branch lengths, no labels');
test_roundtrip('(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;','wkp branch lengths, including root');
test_roundtrip('(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);','wkp distances and leaf names');
test_roundtrip('(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;','wkp distances and all names');
test_roundtrip('((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;','wkp rooted on leaf node');

# From the PHYLIP site:
test_roundtrip('(B,(A,C,E),D);','phylip simple tree');
test_roundtrip('(,(,,),);','phylip no labels');
test_roundtrip('(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);','phylip w/ branch lengths');
test_roundtrip('(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);','phylip w/ internal label');
test_roundtrip('((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700,seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);','phylip raccoon tree');
test_roundtrip('(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;','phylip mammal tree');
test_roundtrip('(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);','phylip mammal tree w/ underbars');
test_roundtrip('A;','phylip single node');
test_roundtrip('((A,B),(C,D));','phylip_quartet');
test_roundtrip('(Alpha,Beta,Gamma,Delta,,Epsilon,,,);','phylip greek');

sub test_roundtrip {
  my $string = shift;
  my $desc = shift;

  my $in = Bio::TreeIO->new(-format => 'newick',
                            -string => $string,
                            -verbose => $verbose
                            );
  my $out = '';
  eval {
    my $t = $in->next_tree;
    $out = $t->as_text('newick');
  };
  return is($out,$string,$desc);
}

sub read_file {
  my $file = shift;
  open my $IN, '<', $file or die "Could not read file '$file': $!\n";
  my (@lines) = <$IN>;
  close $IN;

  @lines = map {$_ =~ s/\\n//g} @lines;
  return join("",@lines);
}
