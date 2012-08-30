#-*-perl-*-
# $Id$

use strict;

use Bio::Root::Test;
use Bio::Tree::Tree;
use Bio::TreeIO;
test_begin( -tests=>125,
	    -requires_modules => [qw(Bio::Phylo)]);

use_ok('Bio::NexmlIO');
diag("WARNING: NeXML parsing for NeXML v0.9 is currently very experimental support");
SKIP: {
    skip("NeXML parsing for NeXML v0.9 is currently very experimental support", 124);
#Read in Data
my $in_nexmlIO = Bio::NexmlIO->new(-file => test_input_file('characters+trees.nexml.xml'));

	#Read in some alignments
	my $aln1 = $in_nexmlIO->next_aln();#, 'nexml matrix to aln' );
	isa_ok($aln1, 'Bio::SimpleAlign', 'obj ok');
	is ($aln1->id,	'DNA sequences', 'aln id');
	my $num =0;
	my @expected_seqs = ('ACGCTCGCATCGCATC', 'ACGCTCGCATCGCATT', 'ACGCTCGCATCGCATG');
	#checking sequence objects
	foreach my $seq_obj ($aln1->each_seq()) {
		$num++;
		
		is( $seq_obj->alphabet, 'dna', "alphabet" );
		is( $seq_obj->display_id, "dna_seq_$num", "display_id");
		is( $seq_obj->seq, $expected_seqs[$num-1], "sequence correct");
	}
	my $aln2 = $in_nexmlIO->next_aln();
	my @alns1;
	push @alns1, $aln1;
	push @alns1, $aln2;
	#checking taxa object
	my %expected_taxa = (dna_seq_1 => 'Homo sapiens', dna_seq_2 => 'Pan paniscus', dna_seq_3 => 'Pan troglodytes');
	my @feats = $aln1->get_all_SeqFeatures();
	foreach my $feat (@feats) {
		if ($feat->has_tag('taxa_id')){
			is ( ($feat->get_tag_values('taxa_id'))[0], 'taxa1', 'taxa id ok' );
			is ( ($feat->get_tag_values('taxa_label'))[0], 'Primary taxa block', 'taxa label ok');
			is ( $feat->get_tag_values('taxon'), 5, 'Number of taxa ok')
		}
		else{
			my $seq_num = ($feat->get_tag_values('id'))[0];
			is ( ($feat->get_tag_values('taxon'))[0], $expected_taxa{$seq_num}, "$seq_num taxon ok" )
		}
	}
	
	#Read in some sequences
	ok( my $seq1 = $in_nexmlIO->next_seq() );
	isa_ok($seq1, 'Bio::Seq');
	is( $seq1->alphabet,		'dna',					"alphabet" );
	is( $seq1->primary_id,	'dna_seq_1',	"primary_id");
	is( $seq1->display_id,	'dna_seq_1',			"display_id");
	is( $seq1->seq,			'ACGCTCGCATCGCATC',		"sequence");

	#checking second sequence object
	ok( my $seq2 = $in_nexmlIO->next_seq() );
	is( $seq2->alphabet,		'dna',					"alphabet" );
	is( $seq2->primary_id,	'dna_seq_2',	"primary_id");
	is( $seq2->display_id,	'dna_seq_2',			"display_id");
	is( $seq2->seq,			'ACGCTCGCATCGCATT',		"sequence");
	ok( my $seq3 = $in_nexmlIO->next_seq() );
	ok( my $seq4 = $in_nexmlIO->next_seq() );
	my @seqs1;
	push @seqs1, $seq1;
	push @seqs1, $seq2;
	push @seqs1, $seq3;
	push @seqs1, $seq4;
	
	#Read in some trees
	ok( my $tree1 = $in_nexmlIO->next_tree() );
	isa_ok($tree1, 'Bio::Tree::Tree');
	is( $tree1->get_root_node()->id(), 'n1', "root node");
	my @nodes = $tree1->get_nodes();
	is( @nodes, 9, "number of nodes");
	ok ( my $node7 = $tree1->find_node('n7') );
	is( $node7->branch_length, 0.3247, "branch length");
	is( $node7->ancestor->id, 'n3');
	is( $node7->ancestor->branch_length, '0.34534');
	#Check leaf nodes and taxa
	my %expected_leaves = (
							'n8'	=>	'bird',
							'n9'	=>	'worm',
							'n5'	=>	'dog',
							'n6'	=>	'mouse',
							'n2'	=>	'human'
	);
	ok( my @leaves = $tree1->get_leaf_nodes() );
	is( @leaves, 5, "number of leaf nodes");
	foreach my $leaf (@leaves) {
		my $leafID = $leaf->id();
		ok( exists $expected_leaves{$leaf->id()}, "$leafID exists"  );
		is( $leaf->get_tag_values('taxon'), $expected_leaves{$leaf->id()}, "$leafID taxon");
	}
	my $tree2 = $in_nexmlIO->next_tree();
	my @trees1;
	push @trees1, $tree1;
	push @trees1, $tree2;


#Write Data
diag('Begin tests for write/read roundtrip');
my $outdata = test_output_file();


my $nexml_out = Bio::NexmlIO->new(-file => ">$outdata", -format => 'Nexml');	

ok( $nexml_out->write(-seqs => \@seqs1, -alns =>\@alns1, -trees => \@trees1), "write to stream" );
close($outdata);

#Read in the out file to test roundtrip
my $in_nexmlIO_roundtrip = Bio::NexmlIO->new(-file => $outdata);

 
	#Read in some alignments
	my $aln3 = $in_nexmlIO_roundtrip->next_aln();#, 'nexml matrix to aln' );
	isa_ok($aln3, 'Bio::SimpleAlign', 'obj ok');
	is ($aln3->id,	'DNA sequences', 'aln id');
	$num =0;
	#checking sequence objects
	foreach my $seq_obj ($aln3->each_seq()) {
		$num++;
		
		is( $seq_obj->alphabet, 'dna', "alphabet" );
		is( $seq_obj->display_id, "dna_seq_$num", "display_id");
		is( $seq_obj->seq, $expected_seqs[$num-1], "sequence correct");
	}
	#checking taxa object
	my @feats_r = $aln3->get_all_SeqFeatures();
	foreach my $feat (@feats_r) {
		if ($feat->has_tag('taxa_id')){
			is ( ($feat->get_tag_values('taxa_id'))[0], 'taxa1', 'taxa id ok' );
			is ( ($feat->get_tag_values('taxa_label'))[0], 'Primary taxa block', 'taxa label ok');
			is ( $feat->get_tag_values('taxon'), 5, 'Number of taxa ok')
		}
		else{
			my $seq_num = ($feat->get_tag_values('id'))[0];
			is ( ($feat->get_tag_values('taxon'))[0], $expected_taxa{$seq_num}, "$seq_num taxon ok" )
		}
	}
	#check extract_alns method
	my $alns_outfile = test_output_file();
	ok ( $in_nexmlIO_roundtrip->extract_alns(-file => ">$alns_outfile", -format => "fasta"), 'extract_alns write' );
	close($alns_outfile);
	my $alnIO = Bio::SeqIO->new(-file => "$alns_outfile", -format => 'fasta');
	my $alns_array = $in_nexmlIO_roundtrip->{_seqs};
	my $alnNum = 1;
	while (my $aln = $alnIO->next_seq()) {
		is( $aln->seq, $alns_array->[$alnNum-1]->seq, "extract_alns roundtrip $alnNum" );
		$alnNum++;
	}
	
	#Read in some sequences
	ok( my $seq5 = $in_nexmlIO_roundtrip->next_seq() );
	isa_ok($seq5, 'Bio::Seq');
	is( $seq5->alphabet,		'dna',					"alphabet" );
	is( $seq5->primary_id,	'dna_seq_1',	"primary_id");
	is( $seq5->display_id,	'dna_seq_1',			"display_id");
	is( $seq5->seq,			'ACGCTCGCATCGCATC',		"sequence");

	#checking second sequence object
	ok( my $seq6 = $in_nexmlIO_roundtrip->next_seq() );
	is( $seq6->alphabet,		'dna',					"alphabet" );
 	is( $seq6->primary_id,	'dna_seq_2',	"primary_id");
	is( $seq6->display_id,	'dna_seq_2',			"display_id");
	is( $seq6->seq,			'ACGCTCGCATCGCATT',		"sequence");
	#check extract_seqs method
	my $seqs_outfile = test_output_file();
	ok ( $in_nexmlIO_roundtrip->extract_seqs(-file => ">$seqs_outfile", -format => "fasta"), 'extract_seqs write' );
	close($seqs_outfile);
	my $seqIO = Bio::SeqIO->new(-file => "$seqs_outfile", -format => 'fasta');
	my $seqs_array = $in_nexmlIO_roundtrip->{_seqs};
	my $seqNum = 1;
	while (my $seq = $seqIO->next_seq()) {
		is( $seq->seq, $seqs_array->[$seqNum-1]->seq, "extract_seqs roundtrip $seqNum" );
		$seqNum++;
	}
	
	#Read in some trees
	ok( my $tree3 = $in_nexmlIO_roundtrip->next_tree() );
	isa_ok($tree3, 'Bio::Tree::Tree');
	is( $tree3->get_root_node()->id(), 'n1', "root node");
	my @nodes3 = $tree3->get_nodes();
	is( @nodes3, 9, "number of nodes");
	ok ( my $node7_r = $tree3->find_node('n7') );
	is( $node7_r->branch_length, 0.3247, "branch length");
	is( $node7_r->ancestor->id, 'n3');
	is( $node7_r->ancestor->branch_length, '0.34534');
	#Check leaf nodes and taxa
	ok( my @leaves3 = $tree3->get_leaf_nodes() );
	is( @leaves3, 5, "number of leaf nodes");
	foreach my $leaf (@leaves3) {
		my $leafID = $leaf->id();
		ok( exists $expected_leaves{$leaf->id()}, "$leafID exists"  );
		is( $leaf->get_tag_values('taxon'), $expected_leaves{$leaf->id()}, "$leafID taxon");
	}
	#check extract_trees method
	my $trees_outfile = test_output_file();
	ok ( $in_nexmlIO_roundtrip->extract_trees(-file => ">$trees_outfile", -format => "nexus"), 'extract_trees write' );
	close($seqs_outfile);
	my $treeIO = Bio::TreeIO->new(-file => "$trees_outfile", -format => 'nexus');
	my $trees_array = $in_nexmlIO_roundtrip->{_trees};
	my $treeNum = 1;
	while (my $tree = $treeIO->next_tree()) {
		is( $tree->id, $trees_array->[$treeNum-1]->id, "extract_trees roundtrip $treeNum" );
		$treeNum++;
	}
}
