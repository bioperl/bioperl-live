# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
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
    plan tests => 13;

}
use Bio::Tree::Node;
use Bio::Tree::PhyloNode;

ok(1);

my $node = new Bio::Tree::Node(-leaf => 0,
			       -parent => undef);

ok(! $node->is_leaf);
ok($node->get_parent, undef);

my $node1 = new Bio::Tree::PhyloNode(-leaf => 1,
				     -parent => $node,
				     -bootstrap => 0.25,
				     -id => 'ADH_BOV',
				     -desc => 'Taxon 1');

ok($node1->is_leaf);
ok($node1->get_parent, $node);
ok($node1->id, 'ADH_BOV');
ok($node1->bootstrap, 0.25);
ok($node1->description, 'Taxon 1');

my $node2 = new Bio::Tree::PhyloNode(-leaf => 1,
				     -parent => $node,
				     -bootstrap => 0.30,
				     -id => 'ADH_MUS',
				     -desc => 'Taxon 2');


ok($node2->is_leaf);
ok($node2->get_parent, $node);
ok($node2->id, 'ADH_MUS');
ok($node2->bootstrap, 0.30);
ok($node2->description, 'Taxon 2');
