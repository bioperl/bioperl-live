use strict;
use warnings;
my %ASSEMBLY_TESTS;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 892,
                -requires_module => 'DB_File' );

    use_ok 'Bio::Seq';
    use_ok 'Bio::LocatableSeq';
    use_ok 'Bio::Seq::Quality';
    use_ok 'Bio::Assembly::IO';
    use_ok 'Bio::Assembly::Singlet';
}
use Bio::Root::IO;

#
# Testing Singlet
#
my ($seq_id, $seq_str, $qual, $start, $end, $strand) = ('seq1', 'CAGT-GGT',
  '0 1 2 3 4 5 6 7', 1, 7, -1, );

my $seq = Bio::PrimarySeq->new(-seq => $seq_str, -id  => $seq_id);
my $singlet_id = 'singlet1';
ok my $singlet = Bio::Assembly::Singlet->new( -id => $singlet_id, -seqref =>
  $seq), 'singlet from Bio::PrimarySeq';
isa_ok $singlet, 'Bio::Assembly::Contig';
isa_ok $singlet, 'Bio::Assembly::Singlet';
is $singlet->id, $singlet_id;
ok my $consensus = $singlet->get_consensus_sequence;
isa_ok $consensus, 'Bio::LocatableSeq';
is $consensus->seq, $seq_str;
is $consensus->start, 1;
is $consensus->end, 7;
is $consensus->strand, 1;
is $singlet->get_consensus_length, 8;
is $singlet->get_consensus_quality, undef;
ok my $seqref = $singlet->seqref;
is $seqref->id, $seq_id;
is $seqref->seq, $seq_str;
is $seqref->start, 1;
is $seqref->end, 7;
is $seqref->length, 8;

$seq = Bio::Seq::Quality->new(-seq => $seq_str, -id => $seq_id, -qual => $qual);
ok $singlet = Bio::Assembly::Singlet->new( -id => $singlet_id, -seqref => $seq),
  'singlet from Bio::Seq::Quality';
is $singlet->get_consensus_length, 8;
ok $consensus = $singlet->get_consensus_quality;
isa_ok $consensus, 'Bio::Seq::QualI';
is join(' ', @{$consensus->qual}), $qual;

$seq = Bio::LocatableSeq->new( -seq => $seq_str, -id => $seq_id, -start =>
  $start, -end => $end, -strand => $strand );
ok $singlet = Bio::Assembly::Singlet->new( -id => $singlet_id, -seqref => $seq),
  'singlet from LocatableSeq';
ok $consensus = $singlet->get_consensus_sequence;
is $consensus->start, 1;
is $consensus->end, 7;
is $consensus->strand, -1;
ok $seqref = $singlet->seqref;
is $seqref->start, 1;
is $seqref->end, 7;
is $seqref->strand, -1;

($start, $end) = (20, 26);
$seq = Bio::LocatableSeq->new( -seq => $seq_str, -id => $seq_id, -start =>
  $start, -end => $end, -strand => $strand );
ok $singlet = Bio::Assembly::Singlet->new( -id => $singlet_id, -seqref => $seq),
  'singlet from LocatableSeq with set coordinates';
ok $consensus = $singlet->get_consensus_sequence;
is $consensus->start, 1;
is $consensus->end, 7;
is $consensus->strand, -1;
ok $seqref = $singlet->seqref;
is $seqref->start, 1;
is $seqref->end, 7;
is $seqref->strand, -1;

#
# Testing Contig
#

#
# Testing IO
#

# ACE variants
ok my $aio = Bio::Assembly::IO->new(
    -file   => test_input_file('assembly_with_singlets.ace'),
    -format => 'ace-consed',
);
is $aio->variant, 'consed', 'consed';
ok $aio = Bio::Assembly::IO->new(
    -file   => test_input_file('assembly_with_singlets.ace'),
    -format => 'ace',
);
is $aio->variant, 'consed';
ok $aio->variant('454');
is $aio->variant, '454';

#
# Some PHRAP input
#

my $in = Bio::Assembly::IO->new(
    -file    => test_input_file('consed_project','edit_dir','test_project.phrap.out'),
    -verbose => -1,
);
isa_ok $in, 'Bio::Assembly::IO';
while (my $contig = $in->next_contig) {
    isa_ok $contig, 'Bio::Assembly::Contig';
}

$in = Bio::Assembly::IO->new(
    -file    => test_input_file('consed_project','edit_dir','test_project.phrap.out'),
    -verbose => -1,
);
isa_ok $in, 'Bio::Assembly::IO';
my $sc;
TODO: {
    local $TODO = "phrap parser doesn't include the sequence string in the sequence objects.";
    $in->verbose(2);
    eval {$sc = $in->next_assembly};
    ok !$@;
}

$in->verbose(-1);
$in = Bio::Assembly::IO->new(
    -file    => test_input_file('consed_project','edit_dir','test_project.phrap.out'),
    -verbose => -1,
);
ok $sc = $in->next_assembly;
isa_ok $sc, 'Bio::Assembly::Scaffold';


#
# Testing Scaffold
#

is $sc->id, "NoName";
is $sc->id('test'), "test";

isa_ok $sc->annotation, 'Bio::AnnotationCollectionI';
is $sc->annotation->get_all_annotation_keys, 0,"no annotations in Annotation collection?";
is $sc->get_nof_contigs, 1;
is $sc->get_nof_sequences_in_contigs, 2;
is $sc->get_nof_singlets, 2, "get_nof_singlets";
is $sc->get_contig_seq_ids, 2, "get_contig_seq_ids";
is $sc->get_contig_ids, 1, "get_contig_ids";
is $sc->get_singlet_ids, 2, "get_singlet_ids";

my @phrap_contigs = $sc->all_contigs();
isa_ok $phrap_contigs[0], "Bio::Assembly::Contig",'the contig is a Bio::Assembly::Contig';
my @singlets = $sc->all_singlets();
isa_ok $singlets[0], "Bio::Assembly::Contig", 'the singlet is a Bio::Assembly::Contig';
isa_ok $singlets[0], "Bio::Assembly::Singlet", 'the singlet is a Bio::Assembly::Singlet';

my @contig_seq_ids;
ok @contig_seq_ids = $sc->get_contig_seq_ids, "get_contig_seq_ids";
is @contig_seq_ids, 2;
for my $contig_seq_id (@contig_seq_ids) {
  ok not $contig_seq_id =~ m/contig/i;
}
my @contig_ids;
ok @contig_ids = $sc->get_contig_ids, "get_contig_ids";
is @contig_ids, 1;
for my $contig_id (@contig_ids) {
  ok not $contig_id =~ m/contig/i;
}
my @singlet_ids;
ok @singlet_ids = $sc->get_singlet_ids, "get_singlet_ids";
is @singlet_ids, 2;
for my $singlet_id (@singlet_ids) {
  ok not $singlet_id =~ m/contig/i;
}
my @all_seq_ids;
ok @all_seq_ids = $sc->get_all_seq_ids, "get_all_seq_ids";
for my $seq_id (@all_seq_ids) {
  ok not $seq_id =~ m/contig/i;
}
is @all_seq_ids, 4;

#
# Testing ContigAnalysis
#

#
# Testing ACE
#

# ACE Consed variant (default)
$aio = Bio::Assembly::IO->new(
    -file   => test_input_file('consed_project','edit_dir','test_project.fasta.screen.ace.2'),
    -format => 'ace',
);

my $assembly = $aio->next_assembly();

my @contigs = $assembly->all_contigs();

my $direction = $contigs[0]->strand;
is $direction, 1;

my $features =  $contigs[0]->get_features_collection;

my @contig_features = $features->features;
is @contig_features, 61, 'contig features'; # 59 contig features + 2 seqfeatures

my @annotations = $features->get_features_by_type('Annotation');
is @annotations, 2;

my $had_tag = 0;
for my $an (@annotations) {
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
is $assembly->get_nof_singlets, 0, "get_nof_singlets";
is $assembly->get_contig_seq_ids, 2, "get_contig_seq_ids";
is $assembly->get_contig_ids, 1, "get_contig_ids";
is $assembly->get_singlet_ids, 0, "get_singlet_ids";

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
my $ace_contig = $ace_contigs[0];
isa_ok $ace_contig, "Bio::Assembly::Contig",'the contig is a Bio::Assembly::Contig';

ok my @test_reads = $ace_contig->get_seq_ids();
is scalar @test_reads, 2;
is $test_reads[0], '5704073';
is $test_reads[1], '5762101';

is $assembly->get_nof_sequences_in_contigs, 6;
is $assembly->get_nof_singlets, 33, "get_nof_singlets";
@singlets = $assembly->all_singlets();
isa_ok $singlets[0], "Bio::Assembly::Contig", 'the singlet is a Bio::Assembly::Contig';
isa_ok $singlets[0], "Bio::Assembly::Singlet", 'the singlet is a Bio::Assembly::Singlet';
ok @contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids";
is @contig_seq_ids, 6;
for my $contig_seq_id (@contig_seq_ids) {
  ok not $contig_seq_id =~ m/contig/i;
}
ok @contig_ids = $assembly->get_contig_ids, "get_contig_ids";
is @contig_ids, 3;
for my $contig_id (@contig_ids) {
  ok $contig_id =~ m/contig/i;
}
ok @singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids";
is @singlet_ids, 33;
for my $singlet_id (@singlet_ids) {
  ok $singlet_id =~ m/contig/i;
}
ok @all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids";
for my $seq_id (@all_seq_ids) {
  ok not $seq_id =~ m/contig/i;
}
is @all_seq_ids, 39;

# bug 2758
ok $aio = Bio::Assembly::IO->new(
    -file   => test_input_file('singlet_w_CT.ace'),
    -format => 'ace',
);

# ACE 454 variant
$aio = Bio::Assembly::IO->new(
    -file   => test_input_file('27-contig_Newbler.ace'),
    -format => 'ace-454',
);
$assembly = $aio->next_assembly();
@contigs = $assembly->all_contigs();
# All read positions should be >0
my $contig = $contigs[0];
my $min_aln_coord = undef;
for my $read ($contig->each_seq) {
   my ($feat) = $contig->get_features_collection->get_features_by_type("_aligned_coord:".$read->id);
   my $aln_coord_start = $feat->location->start;
   if ( (not defined $min_aln_coord) or ($aln_coord_start < $min_aln_coord) ) {
      $min_aln_coord = $aln_coord_start;
   }
}
is $min_aln_coord, 1, '454 ACE variant coordinates check';
# The ends of the consensus should be padded
my $left_pad_length  = 29;
my $sequence_length  = 203;
my $right_pad_length = 81;
my $consensus_length = $left_pad_length + $sequence_length + $right_pad_length;
my $cons_seq  = $contig->get_consensus_sequence->seq;
is length $cons_seq, $consensus_length;
$cons_seq =~ m/^(-*).*?(-*)$/;
is length $1, $left_pad_length, '454 ACE variant consensus check';
is length $2, $right_pad_length;
my $cons_qual = $contig->get_consensus_quality->qual;
is scalar @$cons_qual, $consensus_length;
$cons_qual = join ' ', @{$contig->get_consensus_quality->qual};
my $lpad = $left_pad_length x '0 ';
my $rpad = $right_pad_length x '0 ';
$cons_qual =~ m/^($lpad).*($rpad)$/;
ok defined $1;
ok defined $2;

# Writing ACE files
my $asm_infile  = '27-contig_Newbler.ace';
my $asm_outfile = test_output_file();
my $asm_out = Bio::Assembly::IO->new(
    -file    => ">$asm_outfile",
    -format  =>'ace',
);
my $asm_in;
ok $asm_in = Bio::Assembly::IO->new(
    -file    => test_input_file($asm_infile),
    -format  => 'ace',
    -variant => '454',
)->next_assembly, 'writing in the ACE format';
ok $asm_out->write_assembly( -scaffold => $asm_in, -singlets => 1 );

$asm_infile = 'assembly_with_singlets.ace';
ok $asm_in = Bio::Assembly::IO->new(
    -file   => test_input_file($asm_infile),
    -format => 'ace',
)->next_assembly;
ok $asm_out->write_assembly( -scaffold => $asm_in, -singlets => 1  );

$asm_infile = 'reference_ace.ace';
ok $asm_in = Bio::Assembly::IO->new(
    -file   => test_input_file($asm_infile),
    -format => 'ace',
)->next_assembly;
ok $asm_out->write_assembly( -scaffold => $asm_in, -singlets => 1  );


#
# Testing TIGR format
#

# Importing an assembly

$asm_in = Bio::Assembly::IO->new(
    -file   => test_input_file("sample_dataset.tigr"),
    -format => 'tigr',
);
while (my $obj = $asm_in->next_contig) {
    isa_ok $obj, 'Bio::Assembly::Contig'; # Singlets are contigs too
}

$asm_in = Bio::Assembly::IO->new(
    -file   => test_input_file("sample_dataset.tigr"),
    -format => 'tigr',
);

my $scaf_in = $asm_in->next_assembly;
isa_ok $scaf_in, 'Bio::Assembly::Scaffold';
is $scaf_in->id, 'NoName';
is $scaf_in->get_nof_contigs, 13;
is $scaf_in->get_nof_sequences_in_contigs, 36;
is $scaf_in->get_nof_singlets, 1;
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
is_deeply [sort $scaf_in->get_contig_seq_ids], \@contigseqids;
is_deeply [sort $scaf_in->get_contig_ids],     \@contigids   ;
is_deeply [sort $scaf_in->get_singlet_ids],    \@singletids  ;
isa_ok $scaf_in->get_seq_by_id('sdsu|SDSU1_RFPERU_001_A09.x01.phd.1'),'Bio::LocatableSeq';
$contig = $scaf_in->get_contig_by_id('106');
isa_ok $contig,'Bio::Assembly::Contig';

# check Contig object SeqFeature::Collection
# should add more specific Contig tests...
my @sfs = $contig->get_features_collection->features; # 5 contig features + 2 seqfeatures
is scalar @sfs, 7;
is $sfs[1]->seq_id(), undef; # should this be undef?
ok $contig->get_features_collection->get_features_by_type('_aligned_coord:sdsu|SDSU_RFPERU_006_E04.x01.phd.1');
isa_ok $scaf_in->annotation, 'Bio::AnnotationCollectionI';
is $scaf_in->annotation->get_all_annotation_keys, 0, "no annotations in Annotation collection?";


# Exporting an assembly
$asm_outfile = test_output_file();
$asm_out = Bio::Assembly::IO->new(
    -file   =>  ">$asm_outfile",
    -format => 'tigr',
);
ok $asm_out->write_assembly( -scaffold => $scaf_in), 'writing in the TIGR format';


#
# Testing maq
# /maj
#
my $file = 'test.maq';
ok $aio = Bio::Assembly::IO->new(
    -file   => test_input_file($file),
    -format => 'maq',
), "init maq IO object";
ok $assembly = $aio->next_assembly, "get maq assy";
is $assembly->get_nof_contigs, 11, "got all contigs";
ok open(my $tf, test_input_file($file)), "read test file as text";
my @lines = <$tf>;
is $assembly->get_nof_contig_seqs, scalar @lines, "recorded all maq reads";
ok !$assembly->get_nof_singlets, "no singlets";

ok $aio = Bio::Assembly::IO->new( -file => test_input_file($file),
                                  -format => 'maq' );
isa_ok $aio, 'Bio::Assembly::IO';
while (my $contig = $aio->next_contig) {
    isa_ok $contig, 'Bio::Assembly::Contig';
}

#
# Testing maq with singlets
#
$file = 'test_singlets.maq';
ok $aio = Bio::Assembly::IO->new(
    -file   => test_input_file($file),
    -format => 'maq',
);
ok $assembly = $aio->next_assembly, "get maq assy";
isa_ok $aio, 'Bio::Assembly::IO';


ok @contig_seq_ids = $assembly->get_contig_seq_ids, "get_contig_seq_ids";
is @contig_seq_ids, 246;
for my $contig_seq_id (@contig_seq_ids) {
  ok not $contig_seq_id =~ m/maq_assy/i;
}

ok @contig_ids = $assembly->get_contig_ids, "get_contig_ids";
is @contig_ids, 37;
for my $contig_id (@contig_ids) {
  ok $contig_id =~ m/maq_assy/i;
}

ok @singlet_ids = $assembly->get_singlet_ids, "get_singlet_ids";
is @singlet_ids, 4;
for my $singlet_id (@singlet_ids) {
  ok $singlet_id =~ m/maq_assy/i;
}

ok @all_seq_ids = $assembly->get_all_seq_ids, "get_all_seq_ids";
for my $seq_id (@all_seq_ids) {
  ok not $seq_id =~ m/maq_assy/i;
}
is @all_seq_ids, 250;

ok $aio = Bio::Assembly::IO->new(
    -file   => test_input_file($file),
    -format => 'maq',
);
while (my $contig = $aio->next_contig) {
    isa_ok $contig, 'Bio::Assembly::Contig';
}

##############################################
# test format() and variant() in Bio::RootIO
##############################################

$in = Bio::Assembly::IO->new(
   -file   => test_input_file('assembly_with_singlets.ace'),
);
is $in->format, 'ace';
is $in->variant, 'consed';

exit;
