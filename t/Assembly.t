# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    $error = 0;
    if( $@ ) {
        use lib 't\lib';
    }
    use Test::More;

    $NUMTESTS = 22;
    plan tests => $NUMTESTS;
}

if( $error ==  1 ) {
    exit(0);
}

#syntax test

SKIP: {
    eval { require DB_File };
	skip("DB_File not installed",19) if $@;

	require_ok('Bio::Assembly::IO');
	require_ok('Bio::Assembly::Scaffold');
	require_ok('Bio::Assembly::Contig');
	require_ok('Bio::Assembly::ContigAnalysis');

	use Data::Dumper;
	
	#
	# Testing IO
	#
	
	# -file => ">".Bio::Root::IO->catfile("t","data","primaryseq.embl")
	
	my $in = Bio::Assembly::IO->new
		(-file=>Bio::Root::IO->catfile
		 ("t","data","consed_project","edit_dir","test_project.phrap.out"));
	
	isa_ok($in, 'Bio::Assembly::IO');
	
	my $sc = $in->next_assembly;
	isa_ok($sc, 'Bio::Assembly::Scaffold');

	#print Dumper $sc;
	
	#
	# Testing Scaffold
	#
	
	
	is $sc->id, "NoName";
	is $sc->id('test'), "test";
	
	isa_ok($sc->annotation, 'Bio::AnnotationCollectionI');
	is $sc->annotation->get_all_annotation_keys, 0,"no annotations in Annotation collection?";
	is $sc->get_nof_contigs, 1;
	is $sc->get_nof_sequences_in_contigs, 2;
	
	SKIP : {
		skip("TODO: get_nof_singlets() should return a number", 1);
		is($sc->get_nof_singlets, 1, "get_nof_singlets");
	}
	SKIP : {
		skip("TODO: get_seq_ids() should return a number", 1);
		is($sc->get_seq_ids, 2, "get_seq_ids");
	}
	is($sc->get_contig_ids, 1, "get_contig_ids");
	SKIP: {
		skip("TODO: get_singlet_ids() should return a list", 1);
		is $sc->get_singlet_ids, 0;
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
		-file=>Bio::Root::IO->catfile
		 ("t","data","consed_project","edit_dir","test_project.fasta.screen.ace.2"),
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
}