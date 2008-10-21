# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 47,
	-requires_module => 'DB_File');
	
	use_ok('Bio::Assembly::IO');
}

#
# Testing IO
#


my $in = Bio::Assembly::IO->new
	(-file => test_input_file('consed_project','edit_dir','test_project.phrap.out'));

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
is($sc->get_nof_singlets, 2, "get_nof_singlets");
is($sc->get_contig_seq_ids, 2, "get_contig_seq_ids");
is($sc->get_contig_ids, 1, "get_contig_ids");
is($sc->get_singlet_ids, 2, "get_singlet_ids");

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
	-file=>test_input_file('consed_project','edit_dir','test_project.fasta.screen.ace.2'),
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

is $assembly->get_nof_contigs, 1;
is $assembly->get_nof_sequences_in_contigs, 2;
is($assembly->get_nof_singlets, 0, "get_nof_singlets");
is($assembly->get_contig_seq_ids, 2, "get_contig_seq_ids");
is($assembly->get_contig_ids, 1, "get_contig_ids");
is($assembly->get_singlet_ids, 0, "get_singlet_ids");


$aio = Bio::Assembly::IO->new(
	-file=>test_input_file('assembly_with_singlets.ace'),
	-format=>'ace',
);
$assembly = $aio->next_assembly();
is $assembly->get_nof_contigs, 3;
is $assembly->get_nof_sequences_in_contigs, 6;
is($assembly->get_nof_singlets, 33, "get_nof_singlets");

is($assembly->get_contig_seq_ids, 6, "get_contig_seq_ids");

is($assembly->get_contig_ids, 3, "get_contig_ids");
is($assembly->get_singlet_ids, 33, "get_singlet_ids");




#
# Testing TIGR format
#

# Importing an assembly

my $asm_in = Bio::Assembly::IO->new(
    -file => test_input_file("sample_dataset.tasm "),
    -format=>'tigr'
);
my $scaf_in = $asm_in->next_assembly;

isa_ok($scaf_in, 'Bio::Assembly::Scaffold');
is($scaf_in->id, 'NoName');
is($scaf_in->get_nof_contigs, 13);
is($scaf_in->get_nof_sequences_in_contigs, 36);
is($scaf_in->get_nof_singlets, 0);
my @seqids = sort qw(sdsu|SDSU1_RFPERU_001_A09.x01.phd.1
sdsu|SDSU1_RFPERU_001_B03.x01.phd.1 sdsu|SDSU1_RFPERU_001_B04.x01.phd.1
sdsu|SDSU1_RFPERU_001_E04.x01.phd.1 sdsu|SDSU_RFPERU_002_A01.x01.phd.1
sdsu|SDSU_RFPERU_002_B07.x01.phd.1 sdsu|SDSU_RFPERU_002_C12.x01.phd.1
sdsu|SDSU_RFPERU_002_D08.x01.phd.1 sdsu|SDSU_RFPERU_002_H12.x01.phd.1
sdsu|SDSU_RFPERU_003_G09.x01.phd.1 sdsu|SDSU_RFPERU_004_H12.x01.phd.1
sdsu|SDSU_RFPERU_005_F02.x01.phd.1 sdsu|SDSU_RFPERU_006_D03.x01.phd.1
sdsu|SDSU_RFPERU_006_E04.x01.phd.1 sdsu|SDSU_RFPERU_006_E05.x01.phd.1
sdsu|SDSU_RFPERU_006_H08.x01.phd.1 sdsu|SDSU_RFPERU_007_E09.x01.phd.1
sdsu|SDSU_RFPERU_007_F06.x01.phd.1 sdsu|SDSU_RFPERU_008_B02.x01.phd.1
sdsu|SDSU_RFPERU_009_E07.x01.phd.1 sdsu|SDSU_RFPERU_010_B05.x01.phd.1
sdsu|SDSU_RFPERU_010_B06.x01.phd.1 sdsu|SDSU_RFPERU_010_C09.x01.phd.1
sdsu|SDSU_RFPERU_010_D10.x01.phd.1 sdsu|SDSU_RFPERU_012_H02.x01.phd.1
sdsu|SDSU_RFPERU_013_B05.x01.phd.1 sdsu|SDSU_RFPERU_013_C07.x01.phd.1
sdsu|SDSU_RFPERU_013_C08.x01.phd.1 sdsu|SDSU_RFPERU_013_G10.x01.phd.1
sdsu|SDSU_RFPERU_013_H05.x01.phd.1 sdsu|SDSU_RFPERU_014_H06.x01.phd.1
sdsu|SDSU_RFPERU_015_A05.x01.phd.1 sdsu|SDSU_RFPERU_015_C06.x01.phd.1
sdsu|SDSU_RFPERU_015_E04.x01.phd.1 sdsu|SDSU_RFPERU_015_G04.x01.phd.1
sdsu|SDSU_RFPERU_015_H03.x01.phd.1);
my @contigids = sort qw(106 144 148 17 185 2 210 36 453 500 613 668 93);
is_deeply([sort $scaf_in->get_contig_seq_ids], \@seqids);
is_deeply([sort $scaf_in->get_contig_ids], \@contigids);
is_deeply([$scaf_in->get_singlet_ids], []);
isa_ok($scaf_in->get_seq_by_id('sdsu|SDSU1_RFPERU_001_A09.x01.phd.1'),'Bio::LocatableSeq');
my $contig = $scaf_in->get_contig_by_id('106');
isa_ok($contig,'Bio::Assembly::Contig');
# check Contig object SeqFeature::Collection
# should add more specific Contig tests...
my @sfs = $contig->get_features_collection->get_all_features;
is(scalar(@sfs), 5);
my %primary_tags = map { $_->primary_tag => 1 } @sfs;
ok exists $primary_tags{'_aligned_coord:sdsu|SDSU_RFPERU_006_E04.x01.phd.1'};
is($sfs[1]->seq_id(), undef); # should this be undef?

isa_ok($scaf_in->annotation, 'Bio::AnnotationCollectionI');
is($scaf_in->annotation->get_all_annotation_keys, 0,"no annotations in Annotation collection?");

# Exporting an assembly

my $asm_outfile = test_output_file();
my $asm_out = Bio::Assembly::IO->new(
    -file=> ">$asm_outfile",
    -format=>'tigr'
);

ok $asm_out->write_assembly( -scaffold => $scaf_in);

