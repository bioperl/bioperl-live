# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
  use lib '.';
  use Bio::Root::Test;
  use File::Temp qw(tempfile);
  test_begin( -tests => 40 );
  use_ok('Bio::Tree::Node');
  use_ok('Bio::Tree::AlleleNode');
  use_ok('Bio::TreeIO');
}

my $node1 = Bio::Tree::Node->new();
my $node2 = Bio::Tree::Node->new();
ok( $node1->is_Leaf() );
is( $node1->ancestor, undef );

# tests for tags
ok !$node1->has_tag('test');
is $node1->add_tag_value( 'test', 'a' ), 1;
ok $node1->has_tag('test');
is $node1->add_tag_value( 'test', 'b' ), 2;
my @tags = $node1->get_tag_values('test');
is scalar @tags, 2;
is scalar $node1->get_tag_values('test'), 'a', 'retrieve the first value';

is $node1->remove_tag('test2'), 0;
is $node1->remove_tag('test'),  1;
ok !$node1->has_tag('test');
is $node1->set_tag_value( 'test', ( 'a', 'b', 'c' ) ), 3;
is $node1->remove_all_tags(), undef;
ok !$node1->has_tag('test');

my $pnode = Bio::Tree::Node->new();
$pnode->add_Descendent($node1);
is( $node1->ancestor, $pnode );
$pnode->add_Descendent($node2);
is( $node2->ancestor, $pnode );

ok( !$pnode->is_Leaf );

my $phylo_node = Bio::Tree::Node->new(
  -bootstrap => 0.25,
  -id        => 'ADH_BOV',
  -desc      => 'Taxon 1'
);
$node1->add_Descendent($phylo_node);
ok( !$node1->is_Leaf );
is( $phylo_node->ancestor,    $node1 );
is( $phylo_node->id,          'ADH_BOV' );
is( $phylo_node->bootstrap,   0.25 );
is( $phylo_node->description, 'Taxon 1' );

is $phylo_node->ancestor($node2), $node2;
ok $node1->is_Leaf;
is my @descs = $node2->each_Descendent, 1;
is $descs[0], $phylo_node;

my $allele_node = Bio::Tree::AlleleNode->new();
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm1',
    -alleles     => [0]
  )
);
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm3',
    -alleles     => [ 1, 1 ]
  )
);
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm4',
    -alleles     => [ 0, 4 ]
  )
);
ok($allele_node);
my @mkrs = $allele_node->get_marker_names;

is( @mkrs, 3 );
my ($m3) = $allele_node->get_Genotypes( -marker => 'm3' );
is( $m3->get_Alleles, 2 );
my ($a1) = $allele_node->get_Genotypes( -marker => 'm1' )->get_Alleles;
is( $a1, 0 );

my ( $a2, $a3 ) = $allele_node->get_Genotypes( -marker => 'm4' )->get_Alleles;
is( $a2, 0 );
is( $a3, 4 );

# bug 2877
my $str = "(A:52,(B:46,C:50):11,D:70)68"; 
my $in = Bio::TreeIO->new(
  -internal_node_id => 'bootstrap',
  -format => 'nhx',
  -string => $str,
);
my $t = $in->next_tree;
  my $s;
  my $old_root = $t->get_root_node();
  my ($b) = $t->find_node( -id => "B" );
  my $b_anc = $b->ancestor;

  my $r = $b->create_node_on_branch( -FRACTION => 0.5 );
  $r->id('fake');

  # before reroot
  is( $t->as_text('newick',$in->params), "(A:52,(C:50,(B:23)fake:23):11,D:70)68;", 'with fake node' );

  # after reroot
  $t->reroot($r);
  is( $t->as_text('newick',$in->params), "(B:23,(C:50,(A:52,D:70)68:11):23)fake;",
    "after reroot on fake node" );
  $t->reroot($b);

  is( $t->as_text('newick',$in->params), "(((C:50,(A:52,D:70)68:11):23)fake:23)B;", "reroot on B" );

  $t->reroot($b_anc);
  $t->splice( -remove_id => 'fake' );

  is(
    $t->as_text('newick',$in->params),
    "(B:23,C:50,(A:52,D:70)68:11);",
    "remove fake node, reroot on former B anc"
  );
  $t->reroot($old_root);
  is( $t->as_text('newick',$in->params), "(A:52,(B:23,C:50):11,D:70)68;", "roundtrip" );
