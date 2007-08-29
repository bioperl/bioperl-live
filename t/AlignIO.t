# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 283);
	
	use_ok('Bio::AlignIO');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PSI format  
$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("testaln.psi"),
    '-format'	=> 'psi');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'QUERY/1-798');
is($aln->no_sequences, 56);

# ARP format
$str  = Bio::AlignIO ->new(
    '-file'	=> test_input_file("testaln.arp"),
    '-format'	=> 'arp');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '01/1-399','ARP get_nse()');
is($aln->no_sequences, 60,'ARP no_sequences()');
is($aln->description, 'Mandenka', 'ARP description()');
is($str->datatype, 'DNA', 'ARP SeqIO datatype()');
$str  = Bio::AlignIO->new(
    '-file'	=> test_input_file("testaln2.arp"),
    '-format'	=> 'arp');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '000/1-29','ARP get_nse()');
is($aln->no_sequences, 3,'ARP no_sequences()');
is($aln->description, 'Population 1', 'ARP description()');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->get_nse, '001/1-29','ARP get_nse()');
is($aln->no_sequences, 8,'ARP no_sequences()');
is($aln->description, 'Population 2', 'ARP description()');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->get_nse, '024/1-29','ARP get_nse()');
is($aln->no_sequences, 6,'ARP no_sequences()');
is($aln->description, 'Population 3', 'ARP description()');

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
is($aln->get_seq_by_pos(1)->get_nse, 'L43844.1/2-149');
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
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'AJ295015.1/58-1');
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
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'RAD25_SCHPO/5-240');
is($aln->accession, 'PF00244.9');
is($aln->id, '14-3-3');
is($aln->description,'14-3-3 protein');
($ann) = $aln->annotation->get_Annotations('gathering_threshold');
isa_ok($ann, 'Bio::Annotation::SimpleValue');
is($ann->display_text, '25.00 25.00; 25.00 25.00;', 'Pfam annotation');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, 'COMB_CLOAB/6-235');
is($aln->accession, 'PF04029.4');
is($aln->id, '2-ph_phosp');
is($aln->description,'2-phosphosulpholactate phosphatase');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
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

# PFAM format (no annotation)
$str = Bio::AlignIO->new(
	  '-file' => test_input_file("testaln.pfam"));
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246');

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			    '-format' => 'pfam');
$status = $strout->write_aln($aln);
is($status, 1, " pfam output test");

# MAF
$str = Bio::AlignIO->new(
	  '-file' => test_input_file("humor.maf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'NM_006987/1-5000', "maf input test";
is $aln->get_seq_by_pos(1)->strand, '-';

# MSF
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.msf"));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', "msf input test";

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			    '-format' => 'msf');
$status = $strout->write_aln($aln);
is $status, 1, "msf output test";


# FASTA
$str = Bio::AlignIO->new(
		 -file => test_input_file("testaln.fasta"), 
		 -format => 'fasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', 
  "fasta input test ";
is ($aln->get_seq_by_pos(1)->description, 'DESCRIPTION HERE', 
    "fasta input test for description");
is ($aln->get_seq_by_pos(11)->display_id, 'AK_YEAST',
    "fasta input test for id");

is ($aln->get_seq_by_pos(2)->end, 318,
    "fasta input test for end");

is ($aln->get_seq_by_pos(11)->description, 'A COMMENT FOR YEAST', 
    "fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'fasta');
$status = $strout->write_aln($aln);
is $status, 1,"fasta output test";

my $in = Bio::AlignIO->newFh(
   '-file'  => test_input_file("testaln.fasta"), 
			       '-format' => 'fasta');
my $out = Bio::AlignIO->newFh(
   '-file' => ">".test_output_file(), 
				'-format' => 'pfam');
while ( $aln = <$in>) {
    is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431',
     "filehandle input test  ";
    $status = print $out $aln;
    last;
}
is $status, 1, "filehandle output test";


# SELEX
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.selex"),
			   '-format' => 'selex');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/114-431', "selex format test ";

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'selex');
$status = $strout->write_aln($aln);
is $status, 1, "selex output test";


# MASE
$str = Bio::AlignIO->new(
   '-file' => test_input_file("testaln.mase"),
			   '-format' => 'mase');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'AK1H_ECOLI/1-318', "mase input test ";


# PRODOM
$str = Bio::AlignIO->new(
   '-file' => test_input_file("testaln.prodom"),
			   '-format' => 'prodom');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "prodom input test ";


# CLUSTAL
my $outfile = test_output_file();
$strout = Bio::AlignIO->new(
   '-file' => ">$outfile", 
			      '-format' => 'clustalw');
$status = $strout->write_aln($aln);
is $status, 1, "clustalw (.aln) output test";
undef $strout;
$str = Bio::AlignIO->new(
   '-file'=> $outfile, 
			   '-format' => 'clustalw');
$aln = $str->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'P04777/1-33', "clustalw (.aln) input test";
my $io = Bio::AlignIO->new(
   -file => test_input_file("testaln.aln") );
$aln = $io->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string, "MNEGEHQIKLDELFEKLLRARKIFKNKDVLRHSWEPKDLPHRHEQIEA".
"LAQILVPVLRGETMKIIFCGHHACELGEDRGTKGFVIDELKDVDEDRNGKVDVIEINCEHMDTHYRVLPNIAKLF".
"DDCTGIGVPMHGGPTDEVTAKLKQVIDMKERFVIIVLDEIDKLVKKSGDEVLYSLTRINTELKRAKVSVIGISND".
"LKFKEYLDPRVLSSLSEEEVVFPPYDANQLRDILTQRAEEAFYPGVLDEGVIPLCAALAAREHGDARKALDLLRV".
"AGEIAEREGASKVTEKHVWKAQEKIEQDMMEEVIKTLPLQSKVLLYAIVLLDENGDLPANTGDVYAVYRELCEYI".
"DLEPLTQRRISDLINELDMLGIINAKVVSKGRYGRTKEIRLMVTSYKIRNVLRYDYSIQPLLTISLKSEQRRLI",
"clustalw consensus_string test";

# BL2SEQ
$str = Bio::AlignIO->new(
   '-file'   => test_input_file("bl2seq.out"),
			 '-format' => 'bl2seq',
			 '-report_type' => 'blastp');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(2)->get_nse, 'ALEU_HORVU/60-360', 
    "BLAST bl2seq format test";

# PHYLIP interleaved
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.phylip"),
    '-format' => 'phylip');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new(
    '-file'  => ">".test_output_file(),
    '-format' => 'phylip');
$status = $strout->write_aln($aln);
is $status, 1, "phylip output test";


# METAFASTA
#print STDERR "Better Metafasta tests needed\n" if $DEBUG;
$io = Bio::AlignIO->new(-verbose => $DEBUG ? $DEBUG : -1, 
   -file => test_input_file('testaln.metafasta'));
$aln = $io->next_aln;
isa_ok($aln,'Bio::Align::AlignI');
is $aln->consensus_string,'CDEFHIJKLMNOPQRSTUVWXYZ', "consensus_string on metafasta";
is $aln->symbol_chars,'23',"symbol_chars() using metafasta";

# NEXUS
$str = Bio::AlignIO->new(
   '-file' => test_input_file('testaln.nexus'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapiens/1-45';
$strout = Bio::AlignIO->new('-file'  => ">".test_output_file(),
			    '-format' => 'nexus', );
$status = $strout->write_aln($aln);
is $status, 1, "nexus output test";

$str = Bio::AlignIO->new(
   '-file' => test_input_file('Bird_Ovomucoids.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('basic-ladder.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('Kingdoms_DNA.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file('char-interleave.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file('Primate_mtDNA.nex'),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("char-matrix-spaces.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family4nl.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family7n.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("long-names.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("SPAN_Family8a.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("multiline-intrablock-comment.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("Treebase-chlamy-dna.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("quoted-strings1.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("UnaSmithHIV-both.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("quoted-strings2.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("barns-combined.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("radical-whitespace.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
   '-file' => test_input_file("basic-bush.nex"),
			  '-format' => 'nexus');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
$str = Bio::AlignIO->new(
  '-file' => test_input_file("radical-whitespace_02.nex"),
			  '-format' => 'nexus');


# EMBOSS water
$str = Bio::AlignIO->new('-format' => 'emboss',
		 '-file'   => test_input_file('cysprot.water'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'501.50');
is($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/3-342');
is($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-331');
is(sprintf("%.1f",$aln->overall_percentage_identity),33.8);
is(sprintf("%.1f",$aln->average_percentage_identity),40.1);

is($aln->get_seq_by_pos(1)->start, 3);
is($aln->length,364);


# EMBOSS needle
$str = Bio::AlignIO->new('-format' => 'emboss',
	  '-file'   => test_input_file('cysprot.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'499.50');
is($aln->get_seq_by_pos(1)->get_nse,'PAPA_CARPA/1-345');
is($aln->get_seq_by_pos(2)->get_nse,'CATL_HUMAN/1-333');


# EMBOSS water 2.2.x
$str = Bio::AlignIO->new('-format' => 'emboss',
	 '-file'   => test_input_file('cys1_dicdi.water'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/1-343');
is($aln->get_seq_by_pos(2)->get_nse,'CYS1_DICDI-1/1-343');
is($aln->score,'1841.0');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'CYS1_DICDI/29-343');
is($aln->get_seq_by_pos(2)->get_nse,'ALEU_HORVU/61-360');


# EMBOSS water 2.2.x sparse needle
$str = Bio::AlignIO->new(-verbose => $DEBUG,
	  '-format' => 'emboss',
   	'-file'   => test_input_file('sparsealn.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->score,'18.0');
is(sprintf("%.1f",$aln->overall_percentage_identity), 2.1);
is(sprintf("%.1f",$aln->average_percentage_identity), 38.5);
is($aln->get_seq_by_pos(1)->length, 238);
is($aln->length,238);
is($aln->get_seq_by_pos(1)->get_nse,'KV1K_HUMAN/1-108');
is($aln->get_seq_by_pos(2)->get_nse,'IF1Y_HUMAN/1-143');
is($aln->get_seq_by_pos(1)->seq(), 'DIQMTQSPSTLSVSVGDRVTITCEASQTVLSYLNWYQQK'.
   'PGKAPKLLIYAASSLETGVPSRFSGQGSGTBFTFTISSVZPZBFATYYCQZYLDLPRTFGQGTKVDLKR'.
   '-'x130);
is($aln->get_seq_by_pos(2)->seq(), ('-'x94).'PKNKGKGGK-NRRRGKNENESEKRELVFKE'.
   'DGQEYAQVIKMLGNGRLEALCFDGVKRLCHIRGKLRKKVWINTSDIILVGLRDYQDNKADVILKYNADEAR'.
   'SLKAYGGLPEHAKINETDTFGPGDDDEIQFDDIGDDDEDIDDI');
is($aln->is_flush, 1);


# MEGA
$str = Bio::AlignIO->new('-format' => 'mega',
  	'-file'   => test_input_file("hemoglobinA.meg"));

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(1)->get_nse,'Human/1-141');
is($aln->get_seq_by_pos(2)->get_nse,'Horse/1-144');
$aln->unmatch();
is($aln->get_seq_by_pos(3)->subseq(1,10), 'V-LSAADKGN');

$strout = Bio::AlignIO->new('-format' => 'mega',
	  '-file'   => ">" .test_output_file());

$status = $strout->write_aln($aln);
is $status, 1, "mega output test";


# EMBOSS needle
$str = Bio::AlignIO->new('-format' => 'emboss',
	  '-file'   => test_input_file('gf-s71.needle'));
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'MEDVTLFQFTWRKPI-RLQGEIVYKTSETQTIETNKKDVECVANFQENKEVQTDS-VDNGVGENVKKDITISKEVLNLLYDFVRDDSKVNYDRLLEFHKFDKVALETVQKYHVETRNENIILMISSSSRKTLILFGGISHETFCSHQARALLCSSSTSFSIPLPVCAISAVFYSSTQFILGDVSGNISMCSKDKIIFEKKITDGAVTCLEMCRHGLLSGSDDGNIILWQIGTSGLEKLGGTKLTVSDLSRKIRRSSTSNKPVAIVSMQVYVWPSGEEACVATETGGLYLLTLPTLDYKPLSHQTATSINKILFENQFVAVIYHTSNAAVFNSEGLVDEIPFVATLAVR----------PKLVLF--YTSVCVQDITLNCTSPFREFNNEYNPVIKFSKIRFSADLSVING-FRTSSPNSNN-----------------------------------------------');
is($aln->get_seq_by_pos(1)->seq(), 'MEDVTLHHFRWRKPVENKNGEIVYKTSETQTAEISRKDVECVANFQKSQESQTDDFMQNGVGDGIKKEIRISKEVLGHIYDFLRDDSKVNYDRLLEFHKFDKVSLETVQKYHVETRNENIILMISNSSRKTLILFGGLSHETFCSHQARAVLCSSSTTSSLPLPVCAISAVFYSSTQFLLGDISGNISMWTKEKMIFENKVTDGSVTSLELCRYGLLSGSDDGNVILWKVEESKIEKIEGIKLTVSDLSRKIRRSSTSNKPVAIVSMQV----SGDEVCVATETGGLYLLTLPTLESKPLT-QSATSIFKILYEHPYIAVVYHTSNSAIFNSEGLVDEIPFVATLAVRCGAYFIFSNQSRLIIWSMNTRSTVIDENLNCHS-ICSLSND--------------TLQVLDGDFNLNSQSENSATSESENLRISDLQNLRMLKLQNLRTSEFQNFRTSESQYFKKDNGEL');
is($aln->is_flush(), 1);
is($aln->get_seq_by_pos(1)->get_nse,'gf.s71.44/1-448');
is($aln->get_seq_by_pos(2)->get_nse,'Y50C1A.2/1-406');


# PHYLIP sequential/non-interleaved
$strout = Bio::AlignIO->new('-file'  => test_input_file('noninterleaved.phy'),
			    '-format' => 'phylip');
$aln = $strout->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAA'.
   'AGGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGAATT'.
   'TGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGGTTTATCAAAGTAAGACAGTATGATCAGA'.
   'TACCCATAGAGATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCCACACCTGTCAATATAATTG'.
   'GAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' );


# LARGEMULTIFASTA
$str = Bio::AlignIO->new(
   '-file' => test_input_file('little.largemultifasta'),
                         '-format' => 'largemultifasta');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Human:/1-81', "fasta input test ";
is ($aln->get_seq_by_pos(1)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    "fasta input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'Rat:',
    "fasta input test for id");

is ($aln->get_seq_by_pos(3)->description,
    '72.0:1018606-3386002; 73.0:0-14850845; 74.0:0-83355922; SPECIAL_hsApr2003_3.0:0-414023;',
    "fasta input test for description");

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(),
                            '-format' => 'largemultifasta');
$status = $strout->write_aln($aln);
is $status, 1,"fasta output test";


# POA
# just skip on perl 5.6.0 and earlier as it causes a crash on 
# default perl with OS X 10.2
# fink perl 5.6.0 does not seem to have the problem
# can't figure out what it is so just skip for now
SKIP: {
	skip("skipping due to bug in perl 5.6.0 that comes with OS X 10.2", 10) unless ($^O ne 'darwin' || $] > 5.006);
	
	$str = Bio::AlignIO->new(
			  -file   => test_input_file('testaln.po'),
			  -format => 'po',
			  );
	isa_ok($str, 'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	is $aln->no_sequences, 6;
	
	# output is? i.e. does conversion from clustalw to po give the same alignment?
	$str = Bio::AlignIO->new(
		  '-file'   => test_input_file('testaln.aln'),
		  '-format' => 'clustalw');
	isa_ok($str,'Bio::AlignIO');
	$aln = $str->next_aln();
	isa_ok($aln,'Bio::Align::AlignI');
	$strout = Bio::AlignIO->new(
		 '-file'   => ">" . test_output_file(),
		 '-format' => 'po');
	$status = $strout->write_aln($aln);
	is $status, 1, "po output test";
	
	$str = Bio::AlignIO->new(
		 '-file'   => test_input_file('testaln.po'),
		 '-format' => 'po');
	isa_ok($str,'Bio::AlignIO');
	my $aln2 = $str->next_aln();
	isa_ok($aln2,'Bio::Align::AlignI');
	is $aln2->no_sequences, $aln->no_sequences;
	is $aln2->length, $aln->length;
}


# MEME
# this file has no Strand column
$str = Bio::AlignIO->new(
		-file => test_input_file('test.meme'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,25;
is $aln->no_sequences,4;
is $aln->get_seq_by_pos(3)->seq(),"CCTTAAAATAAAATCCCCACCACCA";
is $aln->get_seq_by_pos(3)->strand,"1";

# this file has a Strand column
$str = Bio::AlignIO->new(
		-file => test_input_file('test.meme2'),
		-format => 'meme');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');is $aln->length,20;
is $aln->no_sequences,8;
is $aln->get_seq_by_pos(8)->seq(),"CCAGTCTCCCCTGAATACCC";
is $aln->get_seq_by_pos(7)->strand,"-1";
is $aln->get_seq_by_pos(6)->strand,"1";

# XMFA
$str = Bio::AlignIO->new(
		 -file => test_input_file("testaln.xmfa"), 
		 -format => 'xmfa');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'chrY/1-598', 
  "xmfa input test ";
is ($aln->get_seq_by_pos(2)->description, undef, 
    "xmfa input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'chr7',
    "xmfa input test for id");
is ($aln->get_seq_by_pos(2)->start, 5000,
    "xmfa input test for end");
is ($aln->get_seq_by_pos(2)->end, 5598,
    "xmfa input test for end");
is ($aln->score, 111, 'xmfa alignment score');

$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'chrY/1000-1060', 
  "xmfa input test ";
is ($aln->get_seq_by_pos(2)->description, undef, 
    "xmfa input test for description");
is ($aln->get_seq_by_pos(3)->display_id, 'chr12',
    "xmfa input test for id");
is ($aln->get_seq_by_pos(2)->start, 6000,
    "xmfa input test for end");
is ($aln->get_seq_by_pos(1)->end, 1060,
    "xmfa input test for end");
is ($aln->score, 11, 'xmfa alignment score');

$strout = Bio::AlignIO->new(
   '-file' => ">".test_output_file(), 
			      '-format' => 'xmfa');
$status = $strout->write_aln($aln);
is $status, 1,"xmfa output test";

my %files = (
	'testaln.phylip' 	=> 'phylip',
	'testaln.psi'	 	=> 'psi',
	'testaln.arp'       => 'arp',
	'rfam_tests.stk'    => 'stockholm',
	'testaln.pfam'      => 'pfam',
	'testaln.msf'       => 'msf',
	'testaln.fasta'     => 'fasta',
	'testaln.selex'     => 'selex',
	'testaln.mase'      => 'mase',
	'testaln.prodom'    => 'prodom',
	'testaln.aln'       => 'clustalw',
	'testaln.metafasta' => 'metafasta',
	'testaln.nexus'     => 'nexus',
	'testaln.po'        => 'po',
	'testaln.xmfa'      => 'xmfa'
 );

# input file handles

while (my ($file, $format) = each %files) {
	my $fhin = Bio::AlignIO->newFh(
	   '-file'  => test_input_file($file), 
					   '-format' => $format);
	my $fhout = Bio::AlignIO->newFh(
	   '-file' => ">".test_output_file(), 
					'-format' => 'clustalw');
	while ( $aln = <$fhin>) {
		cmp_ok($aln->no_sequences, '>=', 2, "input filehandle method test : $format");
        last;
	}
}

# output file handles

for my $format (sort values %files) {
	my $fhin = Bio::AlignIO->newFh(
	   '-file' => test_input_file('testaln.aln'), 
					'-format' => 'clustalw');
	my $fhout = Bio::AlignIO->newFh(
	   '-file'  => '>'.test_output_file(), 
					   '-format' => $format);
	while ( $aln = <$fhin>) {
        $status = print $out $aln;
        last;
    }
    is $status, 1, "filehandle output test : $format";
}