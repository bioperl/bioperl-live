# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin( -tests => 1635,
                -requires_module => 'DB_File' );

    use_ok('Bio::Assembly::IO');
}

#
# Testing IO
#

#
# Some PHRAP input
#

my $in = Bio::Assembly::IO->new
    (-file => test_input_file('consed_project','edit_dir','test_project.phrap.out'));
isa_ok($in, 'Bio::Assembly::IO');
while (my $contig = $in->next_contig) {
    isa_ok($contig, 'Bio::Assembly::Contig');
}

$in = Bio::Assembly::IO->new
    (-file => test_input_file('consed_project','edit_dir','test_project.phrap.out'));
isa_ok($in, 'Bio::Assembly::IO');
my $sc;
TODO: {
    local $TODO = "phrap parser doesn't include the sequence string in the sequence objects.";
    $in->verbose(2);
    eval {$sc = $in->next_assembly};
    ok(!$@);
}

$in->verbose(-1);
$in = Bio::Assembly::IO->new
    (-file => test_input_file('consed_project','edit_dir','test_project.phrap.out'));
ok($sc = $in->next_assembly);
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

my @phrap_contigs = $sc->all_contigs();
isa_ok $phrap_contigs[0], "Bio::Assembly::Contig",'the contig is a Bio::Assembly::Contig';
my @singlets = $sc->all_singlets();
isa_ok $singlets[0], "Bio::Assembly::Contig", 'the singlet is a Bio::Assembly::Contig';
isa_ok $singlets[0], "Bio::Assembly::Singlet", 'the singlet is a Bio::Assembly::Singlet';

my @contig_seq_ids;
ok(@contig_seq_ids = $sc->get_contig_seq_ids, "get_contig_seq_ids");
is(@contig_seq_ids, 2);
for my $contig_seq_id (@contig_seq_ids) {
  ok (not $contig_seq_id =~ m/contig/i);
}
my @contig_ids;
ok(@contig_ids = $sc->get_contig_ids, "get_contig_ids");
is(@contig_ids, 1);
for my $contig_id (@contig_ids) {
  ok (not $contig_id =~ m/contig/i);
}
my @singlet_ids;
ok(@singlet_ids = $sc->get_singlet_ids, "get_singlet_ids");
is(@singlet_ids, 2);
for my $singlet_id (@singlet_ids) {
  ok (not $singlet_id =~ m/contig/i);
}
my @all_seq_ids;
ok(@all_seq_ids = $sc->get_all_seq_ids, "get_all_seq_ids");
for my $seq_id (@all_seq_ids) {
  ok (not $seq_id =~ m/contig/i);
}
is(@all_seq_ids, 4);


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
    -format=>'ace'
);

my $assembly = $aio->next_assembly();

my @contigs = $assembly->all_contigs();

my $direction = $contigs[0]->strand;
is $direction, 1;

my $features =  $contigs[0]->get_features_collection;
my @contig_features = $features->get_all_features;
is @contig_features, 8, 'contig features';

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
    -file   => test_input_file('assembly_with_singlets.ace'),
    -format => 'ace',
);

while (my $obj = $aio->next_contig) {
    isa_ok $obj, 'Bio::Assembly::Contig'; # Singlets are contigs too
}

$aio = Bio::Assembly::IO->new(
    -file   => test_input_file('assembly_with_singlets.ace'),
    -format => 'ace',
);

$assembly = $aio->next_assembly();
is $assembly->get_nof_contigs, 3;
my @ace_contigs = $assembly->all_contigs();
isa_ok $ace_contigs[0], "Bio::Assembly::Contig",'the contig is a Bio::Assembly::Contig';
is $assembly->get_nof_sequences_in_contigs, 6;
is($assembly->get_nof_singlets, 33, "get_nof_singlets");
@singlets = $assembly->all_singlets();
isa_ok $singlets[0], "Bio::Assembly::Contig", 'the singlet is a Bio::Assembly::Contig';
isa_ok $singlets[0], "Bio::Assembly::Singlet", 'the singlet is a Bio::Assembly::Singlet';
ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
is(@contig_seq_ids, 6);
for my $contig_seq_id (@contig_seq_ids) {
  ok (not $contig_seq_id =~ m/contig/i);
}
ok(@contig_ids = $assembly->get_contig_ids, "get_contig_ids");
is(@contig_ids, 3);
for my $contig_id (@contig_ids) {
  ok ($contig_id =~ m/contig/i);
}
ok(@singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids");
is(@singlet_ids, 33);
for my $singlet_id (@singlet_ids) {
  ok ($singlet_id =~ m/contig/i);
}
ok(@all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids");
for my $seq_id (@all_seq_ids) {
  ok (not $seq_id =~ m/contig/i);
}
is(@all_seq_ids, 39);


#
# Testing TIGR format
#

# Importing an assembly

my $asm_in = Bio::Assembly::IO->new(
    -file => test_input_file("sample_dataset.tasm "),
    -format=>'tigr'
);
while (my $obj = $aio->next_contig) {
    isa_ok $obj, 'Bio::Assembly::Contig'; # Singlets are contigs too
}

$asm_in = Bio::Assembly::IO->new(
    -file => test_input_file("sample_dataset.tasm "),
    -format=>'tigr'
);

my $scaf_in = $asm_in->next_assembly;
isa_ok($scaf_in, 'Bio::Assembly::Scaffold');
is($scaf_in->id, 'NoName');
is($scaf_in->get_nof_contigs, 13);
is($scaf_in->get_nof_sequences_in_contigs, 36);
is($scaf_in->get_nof_singlets, 1);
my @contigseqids = sort qw(sdsu|SDSU1_RFPERU_001_A09.x01.phd.1
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
my @contigids     = sort qw(106 144 148 17 185 2 210 36 453 500 613 668 93);
my @singletids    = sort qw(123);
my @singletseqids = sort qw(asdf);
is_deeply([sort $scaf_in->get_contig_seq_ids], \@contigseqids);
is_deeply([sort $scaf_in->get_contig_ids],     \@contigids   );
is_deeply([sort $scaf_in->get_singlet_ids],    \@singletids  );
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
is($scaf_in->annotation->get_all_annotation_keys, 0, "no annotations in Annotation collection?");


# Exporting an assembly

my $asm_outfile = test_output_file();
my $asm_out = Bio::Assembly::IO->new(
    -file=> ">$asm_outfile",
    -format=>'tigr'
);

ok $asm_out->write_assembly( -scaffold => $scaf_in);

#
# Testing maq
# /maj
#
my $file = 'test.maq';
ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -format => 'maq' ), "init maq IO object";
ok $assembly = $aio->next_assembly, "get maq assy";
is( $assembly->get_nof_contigs, 11, "got all contigs");
ok open(my $tf, test_input_file($file)), "read test file as text";
my @lines = <$tf>;
is( $assembly->get_nof_contig_seqs, scalar @lines, "recorded all maq reads");
ok !$assembly->get_nof_singlets, "no singlets";

ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -format => 'maq' );
isa_ok($aio, 'Bio::Assembly::IO');
while (my $contig = $aio->next_contig) {
    isa_ok($contig, 'Bio::Assembly::Contig');
}

#
# Testing maq with singlets
#
$file = 'test_singlets.maq';
ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -format => 'maq' );
ok $assembly = $aio->next_assembly, "get maq assy";
isa_ok($aio, 'Bio::Assembly::IO');

@contig_seq_ids;
ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
is(@contig_seq_ids, 246);
for my $contig_seq_id (@contig_seq_ids) {
  ok (not $contig_seq_id =~ m/maq_assy/i);
}
@contig_ids;
ok(@contig_ids = $assembly->get_contig_ids, "get_contig_ids");
is(@contig_ids, 37);
for my $contig_id (@contig_ids) {
  ok ($contig_id =~ m/maq_assy/i);
}
@singlet_ids;
ok(@singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids");
is(@singlet_ids, 4);
for my $singlet_id (@singlet_ids) {
  ok ($singlet_id =~ m/maq_assy/i);
}
@all_seq_ids;
ok(@all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids");
for my $seq_id (@all_seq_ids) {
  ok (not $seq_id =~ m/maq_assy/i);
}
is(@all_seq_ids, 250);

ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -format => 'maq' );
while (my $contig = $aio->next_contig) {
    isa_ok($contig, 'Bio::Assembly::Contig');
}

SKIP : {

    test_skip(-tests => 828,
	      -requires_module => 'Bio::DB::Sam');

#
# Testing sam
# /maj
#
    my ($aio, $assembly, @contig_seq_ids, @singlet_ids, @contig_ids, @all_seq_ids);
    my $file = 'test.bam';
    my $refdb = 'test.ref.fas';
    ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
				      -refdb => test_input_file($refdb),
				      -format => 'sam' ), "init sam IO object";
    isa_ok($aio, 'Bio::Assembly::IO');
    $aio->_current_refseq_id( ($aio->sam->seq_ids)[0] ); # kludge

    while (my $contig = $aio->next_contig) { 
	isa_ok($contig, 'Bio::Assembly::Contig');
    }
    ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
				      -refdb => test_input_file($refdb),
				      -format => 'sam' ),"reopen";
    ok $assembly = $aio->next_assembly, "get sam assy";
    is( $assembly->get_nof_contigs, 21, "got all contigs"); 
    @contig_seq_ids;

    ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
    is(@contig_seq_ids, 334);
    for my $contig_seq_id (@contig_seq_ids) {
	ok ($contig_seq_id =~ m/^SRR/i);
    }
    @contig_ids;
    ok(@contig_ids = $assembly->get_contig_ids, "get_contig_ids");
    is(@contig_ids, 21);
    for my $contig_id (@contig_ids) {
	ok ($contig_id =~ m/sam_assy/i);
    }
    @singlet_ids;
    ok(@singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids");
    is(@singlet_ids, 35);
    for my $singlet_id (@singlet_ids) {
	ok ($singlet_id =~ m/^SRR/i);
    }
    @all_seq_ids;
    ok(@all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids");
    for my $seq_id (@all_seq_ids) {
	ok ($seq_id =~ m/^SRR/i);
    }
    is(@all_seq_ids, 369);
    
}

SKIP : {

    test_skip(-tests => 755,
	      -requires_modules => qw(Bio::DB::Sam Bio::Tools::Run::Samtools),
	      -requires_executable => 'Bio::Tools::Run::Samtools');

#
# Testing bowtie
#


    my $file = 'test.bowtie';
    my $refdb = 'test.ref.fas';
    ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
				      -index => test_input_file($refdb),
				      -format => 'bowtie' ), "init bowtie IO object";
    isa_ok($aio, 'Bio::Assembly::IO');
    $aio->_current_refseq_id( ($aio->sam->seq_ids)[0] ); # kludge

    while (my $contig = $aio->next_contig) {
	isa_ok($contig, 'Bio::Assembly::Contig');
    }
    ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
				      -index => test_input_file($refdb),
				      -format => 'bowtie' ),"reopen";
    ok $assembly = $aio->next_assembly, "get sam assy";
    is( $assembly->get_nof_contigs, 23, "got all contigs");

    ok(@contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids");
    is(@contig_seq_ids, 312);
    for my $contig_seq_id (@contig_seq_ids) {
	ok ($contig_seq_id =~ m/^SRR/i);
    }
    ok(@contig_ids = $assembly->get_contig_ids, "get_contig_ids");
    is(@contig_ids, 23);
    for my $contig_id (@contig_ids) {
	ok ($contig_id =~ m/sam_assy/i);
    }
    ok(@singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids");
    is(@singlet_ids, 36);
    for my $singlet_id (@singlet_ids) {
	ok ($singlet_id =~ m/^sam_assy/i);
    }
    ok(@all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids");
    for my $seq_id (@all_seq_ids) {
	ok ($seq_id =~ m/^SRR/i);
    }
    is(@all_seq_ids, 348);

}
