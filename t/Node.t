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
    plan tests => 17;
}

use Bio::Tree::Node;
use Bio::Tree::AlleleNode;

ok(1);

my $node1 = new Bio::Tree::Node();
my $node2 = new Bio::Tree::Node();
ok($node1->is_Leaf() );
ok($node1->ancestor, undef);

my $pnode = new Bio::Tree::Node();
$pnode->add_Descendent($node1);
ok($node1->ancestor, $pnode);
$pnode->add_Descendent($node2);
ok($node2->ancestor, $pnode);

ok(! $pnode->is_Leaf);

my $phylo_node = new Bio::Tree::Node(-bootstrap => 0.25,
				     -id => 'ADH_BOV',
				     -desc => 'Taxon 1');
$node1->add_Descendent($phylo_node);
ok(! $node1->is_Leaf);
ok($phylo_node->ancestor, $node1);
ok($phylo_node->id, 'ADH_BOV');
ok($phylo_node->bootstrap, 0.25);
ok($phylo_node->description, 'Taxon 1');

my $allele_node = new Bio::Tree::AlleleNode(-alleles => { 'm1' => [ 0 ],
							  'm2' => [ 1 ],
							  'm3' => [ 0,4] });
ok($allele_node);
my @mkrs = $allele_node->get_marker_names;

ok(@mkrs, 3);
ok($allele_node->get_alleles('m3'), 2);
my ($a1) = $allele_node->get_alleles('m1');
ok($a1, 0);

my ($a2,$a3) = $allele_node->get_alleles('m3');
ok($a2, 0);
ok($a3, 4);
