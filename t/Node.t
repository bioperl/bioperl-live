# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

my $error;
use strict;
BEGIN { 
    $error = 0; 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    plan tests => 22;
	use_ok('Bio::Tree::Node');
	use_ok('Bio::Tree::AlleleNode');
}

my $node1 = new Bio::Tree::Node();
my $node2 = new Bio::Tree::Node();
ok($node1->is_Leaf() );
is($node1->ancestor, undef);

my $pnode = new Bio::Tree::Node();
$pnode->add_Descendent($node1);
is($node1->ancestor, $pnode);
$pnode->add_Descendent($node2);
is($node2->ancestor, $pnode);

ok(! $pnode->is_Leaf);

my $phylo_node = new Bio::Tree::Node(-bootstrap => 0.25,
				     -id => 'ADH_BOV',
				     -desc => 'Taxon 1');
$node1->add_Descendent($phylo_node);
ok(! $node1->is_Leaf);
is($phylo_node->ancestor, $node1);
is($phylo_node->id, 'ADH_BOV');
is($phylo_node->bootstrap, 0.25);
is($phylo_node->description, 'Taxon 1');

is $phylo_node->ancestor($node2), $node2;
ok $node1->is_Leaf;
is my @descs = $node2->each_Descendent, 1;
is $descs[0], $phylo_node;

my $allele_node = new Bio::Tree::AlleleNode();
$allele_node->add_Genotype(new Bio::PopGen::Genotype(-marker_name => 'm1',
						     -alleles=>  [ 0 ]));
$allele_node->add_Genotype(new Bio::PopGen::Genotype(-marker_name => 'm3',
						     -alleles=>  [ 1,1 ]));
$allele_node->add_Genotype(new Bio::PopGen::Genotype(-marker_name => 'm4',
						     -alleles=>  [ 0,4 ]));
ok($allele_node);
my @mkrs = $allele_node->get_marker_names;

is(@mkrs, 3);
my ($m3) = $allele_node->get_Genotypes(-marker => 'm3');
is($m3->get_Alleles, 2);
my ($a1) = $allele_node->get_Genotypes(-marker => 'm1')->get_Alleles;
is($a1, 0);

my ($a2,$a3) = $allele_node->get_Genotypes(-marker => 'm4')->get_Alleles;
is($a2, 0);
is($a3, 4);
