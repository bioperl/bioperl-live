# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    #test_begin(-tests => 37);
    
    use_ok('Bio::DB::Tree::Store');
    use_ok('Bio::TreeIO');
    use_ok('Bio::TreeIO::TreeEventDBBuilder');
}

my $verbose = test_debug();

# create a simple store
my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
				    -create => 1,
                                    -dsn    => 'dbname=test_tree.idx');
isa_ok($dbh, "Bio::DB::Tree::Store");

# create the event handler saving into database
my $handler = Bio::TreeIO::TreeEventDBBuilder->new(-store => $dbh);
# create the TreeIO object
my $treeio = Bio::TreeIO->new(-format  => 'newick',
                              -fh      => \*DATA,
                              -handler => $handler);
ok($treeio,"success creating input stream");
is($treeio->_eventHandler,$handler,"DB event handler installed");

# parse the tree
my $tree = $treeio->next_tree();
ok($tree, "tree parsed");
isa_ok($tree, "Bio::DB::Tree::Tree");
ok($tree->root, "tree has root node");
isa_ok($tree->root,"Bio::DB::Tree::Node");
is($tree->id, "1", "used default id setting");
ok(!$tree->is_rooted, "tree is unrooted");

my @children = $tree->root->children; 
is(scalar(@children), 3, "root has 3 children");
is_deeply([sort {($a||'') cmp ($b||'')} (map { $_->id } @children)],
          [(undef,"Bovine","Rodent")],
          "children have expected ids");
is_deeply([sort {$a <=> $b} (map { $_->branch_length } @children)],
          [(0.54939,0.69395,1.21460)],
          "children have expected branch lengths");
my @nodes = ();
while (@children) {
    my $node = shift(@children);
    push(@children, $node->children);
    push(@nodes,$node);
}
is(scalar(@nodes),11,"tree has 11 non-root nodes");

done_testing();
unlink('test_tree.idx');
exit;

__DATA__
(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147,(P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939,Rodent:1.21460);
