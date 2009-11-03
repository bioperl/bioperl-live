# -*-Perl-*- Test Harness script for Bioperl

use strict;
use warnings;

BEGIN {
	use lib '.';
	use Bio::Root::Test;
  
	test_begin(-tests => 98,
			   -requires_modules => [qw(XML::LibXML XML::LibXML::Reader)]);
	
	use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

SKIP: {
  if (1000*$] < 5008) {
    skip("Reader interface only supported in Perl >= 5.8",96);
  } elsif (XML::LibXML::LIBXML_VERSION() <= 20620) {
    skip("Reader not supported for libxml2 <= 2.6.20",96);
  }
  use_ok('Bio::TreeIO::phyloxml');
  
  if ($verbose) {
	diag("libxml version: ", XML::LibXML::LIBXML_VERSION()); 
  }
  
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
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	is($scientificname->value, 'C. elegans');
	if ($verbose > 0) {
	  diag( "Node C Scientific Name: ",$scientificname->value);
	}
	my ($ac3) = $C->annotation->get_nested_Annotations(-keys=>['scientific_name'], -recursive=>1);
	isa_ok( $ac3, 'Bio::Annotation::Collection');
	($scientificname) = $ac2->get_Annotations('_text');
	is($scientificname->value, 'C. elegans');
	if ($verbose > 0) {
	  diag( "Node C Scientific Name: ",$scientificname->value);
	}
	my ($seq) = @{$C->sequence};
	isa_ok( $seq, 'Bio::SeqI');
	my ($seqac) = $seq->annotation;
	isa_ok( $seqac, 'Bio::Annotation::Collection');
	my ($descac) = $seqac->get_nested_Annotations(-keys=>['desc'], -recursive=>1);
	my ($desc) = $descac->get_Annotations('_text');
	is($desc->value, 'alcohol dehydrogenase');
	if ($verbose > 0) {
	  diag( "Node C Sequence description: ",$desc->value);
	}
	($descac) = $seqac->get_nested_Annotations(-keys=>['desc'], -recursive=>1);
	
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
	ok -s $FILE1;
	if ($verbose > 0) {
	  diag(`cat $FILE1`);
	}
  }
  
  # tree5: homolog relationship and sequence relationship
  # <events> <speciations> <duplications> <symbol> <accession> 
  # <sequence_relation> 
  {
	if ($verbose > 0) {
	  diag("\ntree5: events and relations");
	}
	my $tree = $treeio->next_tree;
	isa_ok($tree, 'Bio::Tree::TreeI');
	if ($verbose > 0) {
	  diag("tree id: ",$tree->id);
	}
	my $node = $tree->get_root_node;
	my ($speciationsac) = $node->annotation->get_nested_Annotations(-keys=>['speciations'], -recursive=>1);
	my ($speciationval) = $speciationsac->get_Annotations('_text');
	is($speciationval->value, '1');
	if ($verbose > 0) {
	  diag("root node speciation event: ", $speciationval->value);
	}
	my @children = ($node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	my @leaves = ();
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
	my ($z) = $leaves[0];
	my $z_seq = $z->sequence->[0];
	isa_ok ($z_seq, 'Bio::SeqI');
	my ($z_id) = $z_seq->annotation->get_nested_Annotations('-keys'=>['id_source'], '-recursive'=>1);
	my ($z_id_text) = $z_id->value;
	my @seq_rels = $z_seq->annotation->get_nested_Annotations('-keys'=>['sequence_relation'], '-recursive'=>1);
	foreach my $rel (@seq_rels) {
	  isa_ok($rel, 'Bio::Annotation::Relation');
	  is ($rel->tagname, 'sequence_relation');
	  my $seqto = $rel->to;
	  isa_ok ($seqto, 'Bio::SeqI');
	  my ($seqto_id) = $seqto->annotation->get_nested_Annotations('-keys'=>['id_source'], '-recursive'=>1);
	  my $seqto_text = $seqto_id->value;
	  if ($verbose > 0) {
		diag( "node ", $z_id_text, " has ", $rel->type, " relation to ", $seqto_text);
	  }
	}
	my ($x) = $leaves[1];
	
  
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my @children = ($tree->get_root_node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	my @leaves = ();
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
	my ($z) = $leaves[0];
	my $z_seq = $z->sequence->[0];
	isa_ok ($z_seq, 'Bio::SeqI');
	my ($z_seqname) = $z_seq->annotation->get_nested_Annotations('-keys'=>['name'], '-recursive'=>1); 
	my ($z_seqname_text) = $z_seqname->get_Annotations('_text');
	is ($z_seqname_text->value, 'NADH-dependent butanol dehydrogenase B');
	my ($z_molseq) = $z_seq->seq;
	is ($z_molseq, 'MVDFEYSIPTRIFFGKDKINVLGRELKKYGSKVLIVYGGGSIKRNGIYDK');
	if ($verbose > 0) {
	  diag("Sequence ", $z_seqname_text->value, " is ", $z_molseq);
	}
	my ($z_seqname_text2) = $treeio->read_annotation('-obj'=>$z_seq, '-path'=>'name');
	is ($z_seqname_text->value, $z_seqname_text2);
	my ($y) = $leaves[1];
	my $y_seq = $y->sequence->[0];
	isa_ok ($y_seq, 'Bio::SeqI');

  # add attribute id_source
  $treeio->add_attribute(
        '-obj' => $z_seq,
        '-attr' => "id_source = \"A\""
        );
  $treeio->add_attribute(
        '-obj' => $y_seq,
        '-attr' => "id_source = \"B\""
        );
 
  # add sequence relation
  $treeio->add_phyloXML_annotation(
          '-obj'=>$tree,
          '-xml'=>'<sequence_relation id_ref_0="A" id_ref_1="B" type="orthology"/>'
          );      
 
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my @children = ($tree->get_root_node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	my @leaves = ();
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
	my ($c) = $leaves[0];
	my ($c_id) = $c->annotation->get_nested_Annotations('-keys'=>['id_source'], '-recursive'=>1);
	my @clade_rels = $c->annotation->get_nested_Annotations('-keys'=>['clade_relation'], '-recursive'=>1);
	foreach my $rel (@clade_rels) {
	  isa_ok($rel, 'Bio::Annotation::Relation');
	  is ($rel->tagname, 'clade_relation');
	  my $nodeto = $rel->to;
	  isa_ok ($nodeto, 'Bio::Tree::AnnotatableNode');
	  my ($nodeto_id) = $nodeto->annotation->get_nested_Annotations('-keys'=>['id_source'], '-recursive'=>1);
	  is ($nodeto_id->value, 'b');
	  my ($nodeto_id2) = $treeio->read_annotation('-obj'=>$nodeto, '-path'=>'id_source', '-attr'=>1);
	  is ($nodeto_id->value, $nodeto_id2);
	  if ($verbose > 0) {
		diag( "node ", $c_id->value, " has ", $rel->type, " relation to ", $nodeto_id->value);
	  }
	}
  
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	diag("property:",$annotations[0]) if $verbose;
	my (@keys) = $annotations[0]->get_all_annotation_keys();
	diag("keys:",@keys) if $verbose;
	my (@value) = $annotations[0]->get_Annotations('_text');
	is($value[0]->value, ' 1200 ');
	if ($verbose > 0) {
	  diag( "Annotation NOAA:depth ",$value[0]->value);
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
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	is($value[0]->value, ' 1200 ');
	if ($verbose > 0) {
	  diag( "Annotation NOAA:depth ",$value[0]->value);
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
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	if ($verbose > 0) {
	  diag("tree id: ",$tree->id);
	}
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
	my ($scientificname) = $A->annotation->get_nested_Annotations('-keys'=>['scientific_name'], '-recursive'=>1);
	my ($scientificname_text) = $scientificname->get_Annotations('_text');
	my ($commonname) = $A->annotation->get_nested_Annotations('-keys'=>['common_name'], '-recursive'=>1); 
	my ($commonname_text) = $commonname->get_Annotations('_text');
	my ($rank) = $A->annotation->get_nested_Annotations('-keys'=>['rank'], '-recursive'=>1); 
	my ($rank_text) = $rank->get_Annotations('_text');
	if ($verbose > 0) {
	  diag("node rank is ", $rank_text->value);
	  diag("node scientific name is ", $scientificname_text->value);
	  diag("node common name is ", $commonname_text->value);
	}
	my ($distribution) = $A->annotation->get_nested_Annotations('-keys'=>['distribution'], '-recursive'=>1);
	my ($desc) = $distribution->get_Annotations('desc');
	my ($desc_text) = $desc->get_Annotations('_text');
	if ($verbose > 0) {
	  diag("node distribution is ", $desc_text->value);
	} 
	
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my $node = $tree->get_root_node;
	my @leaves;
	my @children = ($node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
	my ($D) = $leaves[0];
	my ($point) = $treeio->read_annotation('-obj'=>$D, '-path'=>'distribution/point/geodetic_datum', '-attr'=>1);
	is ($point, 'WGS84');
	my ($lat) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'distribution/point/lat');
	my ($long) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'distribution/point/long');
	my ($alt) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'distribution/point/alt');
	is ($lat, '32.880933');
	is ($long, '-117.217543');
	is ($alt, '104');
	if ($verbose > 0) {
	  diag("node distribution lat: $lat long $long alt $alt");
	}
	
	
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
	my $node = $tree->get_root_node;
	my @leaves;
	my @children = ($node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
  my ($D) = $tree->find_node(-id => 'A');
	my ($dateunit) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/unit', '-attr'=>1);
	my ($datemin) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/minimum' );
	my ($datemax) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/maximum' );
	my ($datevalue) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/value');
	is ($dateunit, 'mya');
	is ($datemin, '416.0');
	is ($datemax, '443.7');
	is ($datevalue, '425');
	if ($verbose > 0) {
	  diag("node date unit: $dateunit min $datemin max $datemax value $datevalue");
	}
  
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
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
  
  # add annotation in phyloxml
	if ($verbose > 0) {
	  diag("test add annotation in phyloXML format");
	}
	my $node = $tree->get_root_node;
	my @leaves;
	my @children = ($node);
	for (@children) {
	  push @children, $_->each_Descendent();
	}
	for (@children) {
	  push @leaves, $_ if $_->is_Leaf;
	}
	my ($D) = $leaves[0];
	isa_ok($D, 'Bio::Tree::AnnotatableNode');
	$treeio->add_phyloXML_annotation(
			  -obj => $D, 
			  -xml => "            <name>D</name>
			  <date unit=\"mya\">
				 <desc>my date</desc>
				 <value>manymany million years</value>
			  </date>
  "
			  );
	my ($dateunit) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/unit', '-attr'=>1);
	my ($datevalue) =  $treeio->read_annotation('-obj'=>$D, '-path'=>'date/value');
	is ($dateunit, 'mya');
	is ($datevalue, 'manymany million years');
  
  # write_tree
	if ($verbose > 0) {
	  diag("\ntest write_tree");
	}
	my $FILE1 = test_output_file();
	my $treeout = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$treeout->write_tree($tree);
	ok -s $FILE1;
	if ($verbose > 0) {
	  diag(`cat $FILE1`);
	}
  }
  
  
  # tree14: convert between nhx and phyloxml
  # convert between nhx-phyloxml
  {
	if ($verbose > 0) {
	  diag("\n test translation between nhx and phyloxml");
	}
	ok my $nhxio = Bio::TreeIO->new(
		-verbose => $verbose,
		-format => 'nhx',
		-file   => test_input_file('test.nhx'));
	my $tree = $nhxio->next_tree;
	isa_ok($tree, 'Bio::Tree::TreeI');
	my $FILE1 = test_output_file();
	my $phyloxmlio = Bio::TreeIO->new(-verbose => $verbose,
		-format => 'phyloxml',
		-file   => ">$FILE1");
	$phyloxmlio->write_tree($tree);
	ok -s $FILE1;
	if ($verbose > 0) {
	  diag(`cat $FILE1`);
	}
  }

  # to-do 1. validation
TODO: {
  local $TODO = 'validation not implemented yet';

  my $xsd = "http://www.phyloxml.org/1.10/phyloxml.xsd";
#  for my $XSD ($xsd, XML::LibXML::Schema->new(location => $xsd)) {
#    {
#      my $reader = new XML::LibXML::Reader(
#  location => "test/schema/demo.xml",
#  Schema => $XSD,
#       );
#      ok($reader->finish, "validate using ".(ref($XSD) ? 'XML::LibXML::Schema' : 'Schema file'));
#    }
#    {
#      my $reader = new XML::LibXML::Reader(
#  location => "test/schema/invaliddemo.xml",
#  Schema => $XSD,
#       );
#      eval { $reader->finish };
#      ok($@, "catch validation error for ".(ref($XSD) ? 'XML::LibXML::Schema' : 'Schema file'));
#    }
#
#  }

} 
  
}
