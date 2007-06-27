# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 19,
			   -requires_module => 'DB_File');
	
	use_ok('Bio::Assembly::IO');
}

#
# Testing IO
#


my $in = Bio::Assembly::IO->new
	(-file => test_input_file("consed_project","edit_dir","test_project.phrap.out"));

isa_ok($in, 'Bio::Assembly::IO');

my $sc = $in->next_assembly;
isa_ok($sc, 'Bio::Assembly::Scaffold');

#
# Testing Scaffold
#


is $sc->id, "NoName";
is $sc->id('test'), "test";

isa_ok($sc->annotation, 'Bio::AnnotationCollectionI');
is $sc->annotation->get_all_annotation_keys, 0,"no annotations in Annotation collection?";
is $sc->get_nof_contigs, 1;
is $sc->get_nof_sequences_in_contigs, 2;

TODO: {
	local $TODO = "get_nof_singlets() should return a number";
	is($sc->get_nof_singlets, 1, "get_nof_singlets");
}
TODO: {
	local $TODO = "get_seq_ids() should return a number";
	is($sc->get_seq_ids, 2, "get_seq_ids");
}
is($sc->get_contig_ids, 1, "get_contig_ids");
TODO: {
	local $TODO = "get_singlet_ids() should return a list";
	isnt $sc->get_singlet_ids, 0;
}
	
#
# Testing Contig
#

#
# Testing ContigAnalysis
#

#
# Testing Ace 
#

my $aio = Bio::Assembly::IO->new(
	-file=>test_input_file("consed_project","edit_dir","test_project.fasta.screen.ace.2"),
	-format=>'ace',
);

my $assembly = $aio->next_assembly();
my @contigs = $assembly->all_contigs();

my $direction = $contigs[0]->strand;
is $direction, 1;

my $features =  $contigs[0]->get_features_collection;
my @contig_features = $features->get_all_features;
is @contig_features, 8;

my @annotations = grep {$_->primary_tag eq 'Annotation'} @contig_features;
is @annotations, 2;
my $had_tag = 0;
foreach my $an (@annotations) {
	if ($an->has_tag('extra_info')) {
		$had_tag++;
		is (($an->get_tag_values('extra_info'))[0], "contig extra\ninfo\n");
	}
	elsif ($an->has_tag('comment')){
		$had_tag++;
		is (($an->get_tag_values('comment'))[0], "contig tag\ncomment\n");
	}
}
is $had_tag, 2;
