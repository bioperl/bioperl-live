# -*-Perl-*- Test Harness script for Bioperl

use strict;
use warnings;

BEGIN {
  use lib 't/lib';
  use BioperlTest;

  test_begin(-tests => 73,
             -requires_modules => [qw(XML::LibXML XML::LibXML::Reader)],
            );
  if (1000*$] < 5008) {
    plan skip_all => "Reader interface only supported in Perl >= 5.8";
    exit;
  } elsif (XML::LibXML::LIBXML_VERSION() <= 20620) {
    plan skip_all => "Reader not supported for libxml2 <= 2.6.20";
    exit;
  }
	use_ok('Bio::TreeIO');
	use_ok('Bio::TreeIO::phyloxml');
  diag("libxml version: ", XML::LibXML::LIBXML_VERSION()); 
}

my $verbose = test_debug();

ok my $treeio = Bio::TreeIO->new(
           -verbose => $verbose,
			     -format => 'phyloxml',
			     -file   => test_input_file('phyloxml_examples.xml'));

# tree1: clade and attribute
# <phylogeny> <clade> <name>
{
  if ($verbose > 0) {
    diag("\ntree1: clade and attribute");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  is($tree->id, 'example from Prof. Joe Felsenstein\'s book "Inferring Phylogenies"');
  is($tree->get_tag_values('description'), 'phyloXML allows to use either a "branch_length" attribute or element to indicate branch lengths.');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
    diag("tree description: ", $tree->get_tag_values('description'));
  }
  is($tree->get_tag_values('rooted'), 'true');
  my @nodes = $tree->get_nodes;
  is(@nodes, 5);
  my ($A) = $tree->find_node('A');
  ok($A);
  is($A->branch_length, '0.102');
  if ($verbose > 0) {
    diag("node A: branch_length ", $A->branch_length);
  }
  is($A->ancestor->id, '');
  is($A->ancestor->branch_length, '0.06');
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree2: branch_length
# <branch_length>
{
  if ($verbose > 0) {
    diag("\ntree2: branch_length");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my @nodes = $tree->get_nodes;
  is(@nodes, 5);
  my $A = $tree->find_node('A');
  ok($A);
  is($A->branch_length, '0.102');
  if ($verbose > 0) {
    diag("node A: branch_length ", $A->branch_length);
  }
  is($A->ancestor->id, '');
  is($A->ancestor->branch_length, '0.06');
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree3: confidence (bootstrap)
# <confidence>
{
  if ($verbose > 0) {
    diag("\ntree3: confidence (bootstrap)");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $AB = $tree->find_node('AB');
  ok($AB);
  is($AB->bootstrap, '89');
  if ($verbose > 0) {
    diag("node AB: bootstrap ", $AB->bootstrap);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree4: species and sequence
# <taxonomy> <scientific_name> <sequence> <annotation>
{
  if ($verbose > 0) {
    diag("\ntree4: taxonomy and sequence");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $C = $tree->find_node('C');
  my ($ac) = $C->annotation->get_Annotations('taxonomy');
  isa_ok( $ac, 'Bio::Annotation::Collection');
  my ($ac2) = $ac->get_Annotations('scientific_name');
  isa_ok( $ac2, 'Bio::Annotation::Collection');
  my ($scientificname) = $ac2->get_Annotations('_text');
  is($scientificname->as_text, 'Value: C. elegans');
  if ($verbose > 0) {
    diag( "Node C Scientific Name: ",$scientificname->as_text);
  }
  my ($ac3) = $C->annotation->get_nested_Annotations(-keys=>['scientific_name'], -recursive=>1);
  isa_ok( $ac3, 'Bio::Annotation::Collection');
  ($scientificname) = $ac2->get_Annotations('_text');
  is($scientificname->as_text, 'Value: C. elegans');
  if ($verbose > 0) {
    diag( "Node C Scientific Name: ",$scientificname->as_text);
  }
  
# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree5: homolog relationship and sequence relationship
# <events> <speciations> <duplications> <symbol> <accession> 
{
  if ($verbose > 0) {
    diag("\ntree5: events");
  }
# <sequence_relation> 
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $node = $tree->get_root_node;
  my @children = ($node);
  for (@children) {
    push @children, $_->each_Descendent();
  }
  my ($A) = $children[0];
  isa_ok($A, 'Bio::Tree::AnnotatableNode');
  my $ac = $A->annotation();
  isa_ok($ac, 'Bio::AnnotationCollectionI');
  if ($verbose > 0) {
    diag($A->to_string());
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree6: detailed sequence data
# <mol_seq> <annotation> <code>
{
  if ($verbose > 0) {
    diag("\ntree6: detailed sequence annotation");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree7: network
# <clade_relation> @id_source & @id_ref
{
  if ($verbose > 0) {
    diag("\ntree7: network using id_source/id_ref");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree8: property elements
# <property>
{
  if ($verbose > 0) {
    diag("\ntree8: property");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my ($A) = $tree->find_node('A');
  isa_ok($A, 'Bio::Tree::AnnotatableNode');
  my ($ac) = $A->annotation();
  isa_ok($ac, 'Bio::AnnotationCollectionI');
  my (@annotations) = $ac->get_Annotations('property');
  isa_ok( $annotations[0], 'Bio::Annotation::Collection');
  diag("property:",$annotations[0]);
  my (@keys) = $annotations[0]->get_all_annotation_keys();
  diag("keys:",@keys);
  my (@value) = $annotations[0]->get_Annotations('_text');
  is($value[0]->as_text, 'Value:  1200 ');
  if ($verbose > 0) {
    diag( "Annotation NOAA:depth stringified ",$value[0]->as_text);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree9: property outside tree topology using id refs
# <property> @id_source @id_ref
{
  if ($verbose > 0) {
    diag("\ntree9: property using id_source/id_ref");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $A = $tree->find_node('A');
  isa_ok($A, 'Bio::Tree::AnnotatableNode');
  my $ac = $A->annotation();
  isa_ok($ac, 'Bio::AnnotationCollectionI');
  my @annotations = $ac->get_Annotations('property');
  isa_ok( $annotations[0], 'Bio::Annotation::Collection');
  my @value = $annotations[0]->get_Annotations('_text');
  is($value[0]->as_text, 'Value:  1200 ');
  if ($verbose > 0) {
    diag( "Annotation NOAA:depth stringified ",$value[0]->as_text);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree10: detailed taxonomy and distribution
# <id> <rank> <uri> <common_name> <distribution>
{
  if ($verbose > 0) {
    diag("\ntree10: taxonomy and distribution");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  my $node = $tree->get_root_node;
  my @leaves;
  my @children = ($node);
  for (@children) {
    push @children, $_->each_Descendent();
  }
  for (@children) {
    push @leaves, $_ if $_->is_Leaf;
  }
  my ($A) = $leaves[0];
  isa_ok($A, 'Bio::Tree::AnnotatableNode');
  my $ac = $A->annotation();
  isa_ok($ac, 'Bio::AnnotationCollectionI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree11: phylogeographic information
# <distribution> <point> <lat> <long> <alt>
{
  if ($verbose > 0) {
    diag("\ntree11: phylogenographic information");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '(((A,B),C),D)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree12: date information
# <date> <desc> <value>
{
  if ($verbose > 0) {
    diag("\ntree12: date");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}

# tree13: alignment outside <phylogeny>
# <align:alignment> <seq> 
{
  if ($verbose > 0) {
    diag("\ntree13: alignment outside  <phylogeny>");
  }
  my $tree = $treeio->next_tree;
  isa_ok($tree, 'Bio::Tree::TreeI');
  if ($verbose > 0) {
    diag("tree id: ",$tree->id);
  }
  my $leaves_string = $tree->simplify_to_leaves_string();
  if ($verbose > 0) {
    diag($leaves_string);
  }
  is($leaves_string, '((A,B),C)');

# write_tree
  if ($verbose > 0) {
    diag("\ntest write_tree");
  }
  my $FILE1 = test_output_file();
  my $treeio = Bio::TreeIO->new(-verbose => $verbose,
      -format => 'phyloxml',
      -file   => ">$FILE1");
  $treeio->write_tree($tree);
  ok -s $FILE1;
  if ($verbose > 0) {
    diag(`cat $FILE1`);
  }
}


