# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 87);

	use_ok('Bio::AlignIO::stockholm');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# STOCKHOLM (multiple concatenated files)
# Rfam
$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("rfam_tests.stk"),
    '-format'	=> 'stockholm');
$strout = Bio::AlignIO->new('-file'  => ">".test_output_file(),
		'-format' => 'stockholm', );

isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->source, 'stockholm');
is($aln->get_seq_by_pos(1)->get_nse, 'Z11765.1/1-89');
is($aln->accession, 'RF00006');
is($aln->id, 'Vault');
is($aln->description,'Vault RNA');
# annotation
my ($ann) = $aln->annotation->get_Annotations('alignment_comment');
isa_ok($ann, 'Bio::Annotation::Comment');
is($ann->text,'This family of RNAs are found as part of the enigmatic vault'.
   ' ribonucleoprotein complex. The complex consists of a major vault protein'.
   ' (MVP), two minor vault proteins (VPARP and TEP1), and several small '.
   'untranslated RNA molecules. It has been suggested that the vault complex '.
   'is involved in drug resistance. We have identified a putative novel vault '.
   'RNA on chromosome 5 EMBL:AC005219.','Stockholm annotation');
is($ann->tagname,'alignment_comment','Stockholm annotation');

# test output
$status = $strout->write_aln($aln);
is $status, 1, "stockholm output test";

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->source, 'stockholm');
is($aln->get_seq_by_pos(1)->get_nse, 'L43844.1/2-149');
is($aln->get_seq_by_pos(1)->version, '1');
is($aln->accession, 'RF00007');
is($aln->id, 'U12');
is($aln->description,'U12 minor spliceosomal RNA');
my @anns = $aln->annotation->get_Annotations('reference');
$ann = shift @anns;
isa_ok($ann, 'Bio::Annotation::Reference', 'Stockholm annotation');
$ann = shift @anns;
is($ann->pubmed,'9149533', 'Stockholm annotation');
is($ann->title,
   'Pre-mRNA splicing: the discovery of a new spliceosome doubles the challenge.',
   'Stockholm annotation');
is($ann->authors,'Tarn WY, Steitz JA;', 'Stockholm annotation');
is($ann->location,'Trends Biochem Sci 1997;22:132-137.', 'Stockholm annotation');
# alignment meta data
my $meta = $aln->consensus_meta;
isa_ok($meta, 'Bio::Seq::MetaI');
my ($name) = $meta->meta_names;
is($name,'SS_cons', 'Rfam meta data');
my $meta_str = $meta->named_meta($name);
is($meta_str, '...<<<<<..........>>>>>........<<<<......<<<<......>>>>>>>>'.
   '<<<<<.......>>>>>...........<<<<<<<...<<<<<<<.....>>>>>>>.>>>>>>>..<<<'.
   '<<<<<<.........>>>>>>>>>...', 'Rfam meta data');
$aln = $str->next_aln();
is($aln->source, 'stockholm');
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'AJ295015.1/58-1');
is($aln->get_seq_by_pos(1)->strand, -1);
is($aln->accession, 'RF00008');
is($aln->id, 'Hammerhead_3');
is($aln->description,'Hammerhead ribozyme (type III)');
# alignment meta data
$meta = $aln->consensus_meta;
isa_ok($meta, 'Bio::Seq::MetaI');
($name) = $meta->meta_names;
is($name,'SS_cons', 'Rfam meta data');
$meta_str = $meta->named_meta($name);
is($meta_str, '.<<<<<<..<<<<<.........>>>>>.......<<<<.....................'.
   '...........>>>>...>>>>>>.', 'Rfam meta data');

# STOCKHOLM (Pfam)
$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("pfam_tests.stk"),
    '-format'	=> 'stockholm');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
is($aln->source, 'stockholm');
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'RAD25_SCHPO/5-240');
is($aln->accession, 'PF00244.9');
is($aln->id, '14-3-3');
is($aln->description,'14-3-3 protein');
($ann) = $aln->annotation->get_Annotations('gathering_threshold');
isa_ok($ann, 'Bio::Annotation::SimpleValue');
is($ann->display_text, '25.00 25.00; 25.00 25.00;', 'Pfam annotation');
$aln = $str->next_aln();
is($aln->source, 'stockholm');
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'COMB_CLOAB/6-235');
is($aln->accession, 'PF04029.4');
is($aln->id, '2-ph_phosp');
is($aln->description,'2-phosphosulpholactate phosphatase');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->source, 'stockholm');
is($aln->get_seq_by_pos(1)->get_nse, 'Y278_HAEIN/174-219');
is($aln->accession, 'PF03475.3');
is($aln->id, '3-alpha');
is($aln->description,'3-alpha domain');
# alignment meta data
$meta = $aln->consensus_meta;
isa_ok($meta, 'Bio::Seq::MetaI');
my %test_data = ('SA_cons'  => '6000320010013274....3372052026033.108303630350385563',
                 'SS_cons'  => 'SCBHHHHHHHHHTSCC....CHHHHHHHHTSTT.CCHHHHHHHHHHHHHSSC',
                 'seq_cons' => 'plTVtclsclhasc......stphLcphLshss.Lupsa+cohpK+lspshs',);
for my $name ($meta->meta_names) {
    ok(exists $test_data{$name}, 'Pfam aln meta data');
    $meta_str = $meta->named_meta($name);
    is($meta_str, $test_data{$name}, 'Pfam aln meta data');
}
%test_data = ();
# sequence meta data
%test_data = ('SA'  => '6000320010013274....3372052026033.108303630350385563',
              'SS'  => 'SCBHHHHHHHHHTSCC....CHHHHHHHHTSTT.CCHHHHHHHHHHHHHSSC');
for my $seq ($aln->each_seq) {
    for my $name ($seq->meta_names) {
        ok(exists $test_data{$name}, 'Pfam seq meta data');
        is($seq->named_meta($name), $test_data{$name}, 'Pfam seq meta data');
    }
}

# sequence-specific alignments
# these are generally DBLinks. However, simple DBLinks are not RangeI, so these
# are now Target (which now is-a DBLink). Since these are per-seq, these are
# stored in a full-length SeqFeature as annotation. For now only sequences which
# have this annotation have a SeqFeature.

my @feats = $aln->get_SeqFeatures;
is(scalar(@feats),6);
isa_ok($feats[0], 'Bio::SeqFeatureI');
isa_ok($feats[0]->entire_seq, 'Bio::Seq::Meta');

my ($link)  = $feats[0]->annotation->get_Annotations('dblink');
isa_ok($link, 'Bio::AnnotationI');
isa_ok($link, 'Bio::Annotation::DBLink');
is($link->database, 'PDB');
is($link->start, '178');
is($link->end, '224');
is($link->primary_id, '1o65');
is($link->target_id, '1o65');  # if not set, defaults to primary_id
is($link->optional_id, 'A');
($link)  = $feats[3]->annotation->get_Annotations('dblink');
is($link->database, 'PDB');
is($link->start, '178');
is($link->end, '224');
is($link->target_id, '1o67');

# bug #3420

my $in = Bio::AlignIO->new(-file => test_input_file('tiny.stk'),
                            -format => 'stockholm');
$aln = $in->next_aln;
is($aln->id, 'NoName');
is($aln->get_seq_by_id('a')->display_id, 'a');
is($aln->get_seq_by_id('b')->display_id, 'b');
