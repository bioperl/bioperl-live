# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 25);
    
    use_ok('Bio::DB::Tree::Tree');
    use_ok('Bio::DB::Tree::Node');
    use_ok('Bio::DB::Tree::Store');
    use_ok('Bio::TreeIO');
}
# create a simple store
my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
				    -create => 1,
                                    -dsn    => 'dbname=test_tree_t.idx');
isa_ok($dbh, "Bio::DB::Tree::Store");
my $in = Bio::TreeIO->new(-format => 'newick',
			  -fh     => \*DATA);
my $t = $in->next_tree;
isa_ok($t, "Bio::Tree::TreeI");
my @nodes = 0;

for my $node ( $t->get_nodes(-order => 'breadth') ) {
#    print $node->id || 
#	join(",",map { $_->id } grep { $_->is_Leaf } $node->get_all_Descendents()), "\n";
    my $parent_id = 0;
    if( $node->ancestor ) {
	$parent_id = $node->ancestor->internal_id;
    }
    my $pk = $dbh->insert_node({'-id'     => $node->id,
				'-parent' => $nodes[$parent_id]});
    ok($pk);
    $nodes[$node->internal_id] = $pk;
}
my $tpk = $dbh->insert_tree({-id => 'tree1',-root=>$nodes[$t->get_root_node->internal_id]});
is($tpk,1);

done_testing();
unlink('test_tree_t.idx');
exit;

__DATA__
((A,(Z,X,Y)),(B,((C,D),(E,(F,G)))))
