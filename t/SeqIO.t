# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
use vars qw($DEBUG $TESTCOUNT $NODOM $NOSAX);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	$TESTCOUNT = 367;
	# interpro uses XML::DOM
	eval {require XML::DOM::XPath};
	if ($@) {
		$TESTCOUNT -= 8;
		$NODOM = 1;
		print STDERR "XML::DOM::XPath not found - skipping interpro tests\n";
	}
	# BSML_SAX uses XML::SAX
	eval {require XML::SAX;
			require XML::SAX::Writer;
			require XML::SAX::Base;
			1;
      };
	if( $@ ) {
		$TESTCOUNT -= 14;
		$NOSAX = 1;
		warn "XML::SAX::Base or XML::SAX or XML::SAX::Writer not found - skipping BSML_SAX tests\n";
	}
	plan tests => $TESTCOUNT;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use Bio::Root::IO;
use Bio::Annotation::Collection;

ok(1);

# Set to -1 for release version, so warnings aren't printed
my $verbosity = $DEBUG ? 1 : -1;

#
# Basic read and/or write tests
#


sub read_write {
	my $format = shift;
	my $seq;

	my $str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","test.$format"),
									  '-format' => $format);
	ok $seq = $str->next_seq();
	print "Sequence 1 of 2 from $format stream:\n", $seq->seq, "\n\n" if  $DEBUG;
	unless ($format eq 'raw') {
		ok $seq->id, 'roa1_drome',"ID for format $format";
		ok $seq->length, 358;
	}

	unless ($format eq 'gcg') { # GCG file can contain only one sequence
		ok $seq = $str->next_seq();
		print "Sequence 2 of 2 from $format stream:\n", $seq->seq, $seq->seq, "\n" if $DEBUG;
	}

	my $out = Bio::SeqIO->new('-file'=> ">". Bio::Root::IO->catfile("t","data","$format.out"),
									  '-format' => $format);
	ok $out->write_seq($seq);
	if ($format eq 'fasta') {
		my $id_type;
		ok($id_type = $out->preferred_id_type('accession.version'), 'accession.version');
	}
}

my @formats = qw(gcg fasta raw pir tab);

foreach my $format (@formats) {
	print "======== $format ========\n" if $DEBUG;
	read_write($format);
}

END {
	map { unlink Bio::Root::IO->catfile("t","data","$_.out") } @formats
}

my ($str, $seq,$ast,$temp,$mf,$ent,$out); # predeclare variables for strict

# PIR testing
$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
                                                         "seqfile.pir"),
                       '-format' => 'pir');
ok $str;
$out = new Bio::SeqIO(-format => 'pir', -fh => \*STDOUT);

while( $seq = $str->next_seq()) {
#    ok($seq->id, qr /^[PF]1/ );
    ok($seq->length > 1);
    $out->write_seq($seq) if $verbosity > 0;
}
$out = undef;

# Genbank
$ast = Bio::SeqIO->new( '-format' => 'GenBank' ,
                        '-file' => Bio::Root::IO->catfile("t","data",
                                                          "roa1.genbank"));
$ast->verbose($verbosity);
my $as = $ast->next_seq();
ok $as->molecule, 'mRNA';
ok $as->alphabet, 'dna';
ok($as->primary_id, 3598416);
my @class = $as->species->classification;
ok $class[$#class],'Eukaryota';

$ast = Bio::SeqIO->new( '-format' => 'genbank' ,
                        '-file' => Bio::Root::IO->catfile("t","data",
                                                          "NT_021877.gbk"));
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok $as->molecule, 'DNA';
ok $as->alphabet, 'dna';
ok($as->primary_id, 37539616);
ok($as->accession_number, 'NT_021877');

my ($cds) = grep { $_->primary_tag eq 'CDS' } $as->get_SeqFeatures();
ok(($cds->get_tag_values('transl_except'))[1],
   '(pos:complement(4224..4226),aa:OTHER)');

# test for a DBSOURCE line
$ast = Bio::SeqIO->new(-format => 'genbank' ,
                       -file => Bio::Root::IO->catfile("t","data",
																		 "BAB68554.gb"));
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok $as->molecule, 'linear';
ok $as->alphabet, 'protein';
ok $as->species->common_name, 'Aldabra giant tortoise';
ok($as->primary_id, 15824047);
my $ac = $as->annotation;
ok defined $ac;
my @dblinks = $ac->get_Annotations('dblink');
ok(scalar @dblinks,1);

# test for multi-line SOURCE
$ast = Bio::SeqIO->new(-format => 'genbank' ,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "NC_006346.gb"));
$as = $ast->next_seq;
ok $as->species->binomial,'Bolitoglossa n. sp.';
@class = $as->species->classification;
ok($class[$#class],'Eukaryota');
ok($as->species->common_name,'mitochondrion Bolitoglossa n. sp. RLM-2004 (mushroomtongue salamander)');

# test for unusual common name
$ast = Bio::SeqIO->new(-format => 'genbank' ,
                       -file => Bio::Root::IO->catfile("t","data",
                                                       "AB077698.gb"));
$as = $ast->next_seq;
ok $as->species->common_name,'Homo sapiens cDNA to mRNA';

# embl
$ast = Bio::SeqIO->new( '-format' => 'embl' ,
                        '-file' => Bio::Root::IO->catfile("t","data",
                                                          "roa1.dat"));
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok defined $as->seq;
ok($as->display_id, 'HSHNCPA1');
ok($as->accession_number, 'X79536');
ok($as->seq_version, 1);
ok($as->version, 1);
ok($as->desc, 'H.sapiens mRNA for hnRNPcore protein A1');
ok($as->molecule, 'RNA');
ok($as->alphabet, 'rna');
ok(scalar $as->all_SeqFeatures(), 4);
ok($as->length, 1198);
ok($as->species->binomial(), 'Homo sapiens');

$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data","test.embl"),
'-format' => 'embl');

$seq = $ent->next_seq();

ok(defined $seq->seq(), 1,
   'failure to read Embl with ^ location and badly split double quotes');
ok(scalar $seq->annotation->get_Annotations('reference'), 3);

$out = Bio::SeqIO->new('-file'=> ">embl.out",
       '-format' => 'embl');

ok($out->write_seq($seq),1,
   'failure to write Embl format with ^ < and > locations');
unlink("embl.out");

# embl with no FT
$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data","test.embl"),
'-format' => 'embl');

$seq = $ent->next_seq();
ok($seq);
ok(lc($seq->subseq(1,10)),'gatcagtaga');
ok($seq->length);

# embl repbase
$ent = Bio::SeqIO->new(-file => Bio::Root::IO->catfile("t","data","BEL16-LTR_AG.embl"), -format => 'embl');
$seq = $ent->next_seq;
ok($seq->display_id,'BEL16-LTR_AG');

# kegg
my $kegg = Bio::SeqIO->new( '-format' => 'kegg' ,
    '-file' => Bio::Root::IO->catfile("t","data","AHCYL1.kegg"));

ok($kegg);
$kegg = $kegg->next_seq();
ok($kegg);
ok($kegg->accession, '10768');
ok($kegg->display_id, 'AHCYL1');
ok($kegg->alphabet, 'dna');
ok($kegg->seq);
ok($kegg->translate->seq);
ok(($kegg->annotation->get_Annotations('description'))[0]->text,
   'S-adenosylhomocysteine hydrolase-like 1 [EC:3.3.1.1]');
ok(($kegg->annotation->get_Annotations('pathway'))[0]->text,
   'Metabolism; Amino Acid Metabolism; Methionine metabolism');

ok( (grep {$_->database eq 'KO'} 
     $kegg->annotation->get_Annotations('dblink'))[0]->comment, 
    'adenosylhomocysteinase' );

ok( (grep {$_->database eq 'PATH'} 
     $kegg->annotation->get_Annotations('dblink'))[0]->primary_id,
    'hsa00271' );

# multifile
$mf = Bio::SeqIO::MultiFile->new( '-format' => 'Fasta' ,
  '-files' =>
  [ Bio::Root::IO->catfile("t","data","multi_1.fa"),
    Bio::Root::IO->catfile("t","data","multi_2.fa")]);

ok defined $mf;
my $count = 0;
eval {
    while( $seq = $mf->next_seq() ) {
$count++;
$temp = $seq->display_id;
    }
};
ok( $count ,12);
$temp = undef;

# swissprot
$ast = Bio::SeqIO->new( '-verbosity' => $verbosity,
                        '-format' => 'swiss' ,
                        '-file' => Bio::Root::IO->catfile("t","data",
                                                          "roa1.swiss"));
$as = $ast->next_seq();

ok defined $as->seq;
ok($as->id, 'ROA1_HUMAN', "id is ".$as->id);
ok($as->primary_id =~ /^Bio::(Primary)?Seq/);
ok($as->length, 371);
ok($as->alphabet, 'protein');
ok($as->division, 'HUMAN');
ok(scalar $as->all_SeqFeatures(), 16);

ok(scalar $as->annotation->get_Annotations('reference'), 11);

($ent, $seq, $out,$as) = undef;

#interpro
unless ($NODOM) {
  my $t_file = Bio::Root::IO->catfile("t","data","test.interpro");
  my $a_in = Bio::SeqIO->new( -FILE => $t_file, -FORMAT => 'interpro'); # ici

  my $seq = $a_in->next_seq();
  ok($seq);
  ok($seq->isa('Bio::Seq::RichSeq'));
  ok(scalar( $seq->get_SeqFeatures() ) == 6);

  my($feat) = $seq->get_SeqFeatures();
  ok($feat->isa('Bio::SeqFeature::Generic'));

  ok($feat->display_name eq 'Retinoblastoma-associated protein, B-box');

  ok($seq = $a_in->next_seq());
  ok(scalar( $seq->get_SeqFeatures() ) == 40);

  ok(!($seq = $a_in->next_seq()));
}

#ace
{
    my $t_file = Bio::Root::IO->catfile("t","data","test.ace");
    my( $before );
    {
        local $/ = undef;
        local *BEFORE;
        open BEFORE, $t_file;
        $before = <BEFORE>;
        close BEFORE;
    }

    my $a_in = Bio::SeqIO->new( -FILE => $t_file, -FORMAT => 'ace');
    my( @a_seq );
    while (my $a = $a_in->next_seq) {
        push(@a_seq, $a);
    }

    ok @a_seq, 3, 'wrong number of sequence objects';

    my $esc_name = $a_seq[1]->display_id;
    ok( $esc_name , 'Name; 4% strewn with \ various / escaped characters',
        "bad unescaping of characters, $esc_name");

    ok $a_seq[0]->alphabet, 'protein', 'alphabets incorrectly detected';
    ok $a_seq[1]->alphabet, 'dna', 'alphabets incorrectly detected';

    my $o_file = Bio::Root::IO->catfile("t","data","test.out.ace");
    my $a_out = Bio::SeqIO->new( -FILE => "> $o_file", -FORMAT => 'ace');
    my $a_out_ok = 1;
    foreach my $a (@a_seq) {
        $a_out->write_seq($a) or $a_out_ok = 0;
    }
    undef($a_out);  # Flush to disk
    ok $a_out_ok,1,'error writing sequence';

    my( $after );
    {
        local $/ = undef;
        local *AFTER;
        open AFTER, $o_file;
        $after = <AFTER>;
        close AFTER;
    }
    unlink($o_file);

    ok( ($before and $after and ($before eq $after)),1,
'test output file differs from input');
}

#
# streaming & Bio::RichSeq creation
#

my $stream = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data",
                                                               "test.genbank"),
                             '-format' => 'GenBank');
$stream->verbose($verbosity);
my $seqnum = 0;
my $species;
my @cl;
my $lasts;
my @ids = qw(DDU63596 DDU63595 HUMBDNF);
my @tids = (44689, 44689, 9606);
my @tnames = ("Dictyostelium discoideum","Dictyostelium discoideum","Homo sapiens");
while($seq = $stream->next_seq()) {
    if($seqnum < 3) {
ok $seq->display_id(), $ids[$seqnum];
$species = $seq->species();
@cl = $species->classification();
ok( $species->binomial(), $tnames[$seqnum],
    'species parsing incorrect for genbank');
ok( $cl[3] ne $species->genus(), 1,
    'genus duplicated in genbank parsing');
ok( $species->ncbi_taxid, $tids[$seqnum] );
    }
    $seqnum++;
    $lasts = $seq;
}
ok $lasts->display_id(), "HUMBETGLOA";
my ($ref) = $lasts->annotation->get_Annotations('reference');
ok($ref->medline, 94173918);
$stream->close();

$stream = Bio::SeqIO->new(
        '-file' => Bio::Root::IO->catfile("t","data","test.genbank.noseq"),
        '-format' => 'GenBank',
        );
$seqnum = 0;
while($seq = $stream->next_seq()) {
    if($seqnum < 3) {
        ok $seq->display_id(), $ids[$seqnum];
    } elsif( $seq->display_id eq 'M37762') {
        ok( ($seq->get_keywords())[0], 'neurotrophic factor');
    }
    $seqnum++;
}
ok $seqnum, 5, "Total number of sequences in test file";

$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data",
                                                          "test.embl"),
                        '-format' => 'embl');
$ent->verbose($verbosity);
$seq = $ent->next_seq();
$species = $seq->species();
@cl = $species->classification();
ok( $cl[3] ne $species->genus(), 1, 'genus duplicated in EMBL parsing');
$ent->close();

$seq = Bio::SeqIO->new( '-format' => 'GenBank' ,
                        '-file' =>Bio::Root::IO->catfile("t","data",
                                                         "testfuzzy.genbank"));
$seq->verbose($verbosity);
ok(defined($as = $seq->next_seq()));

my @features = $as->all_SeqFeatures();
ok(@features,21);
my $lastfeature = pop @features;

# this is a split location; the root doesn't have strand
ok($lastfeature->strand, undef);
my $location = $lastfeature->location;
$location->verbose(-1); # silence the warning of undef seq_id()
# see above; splitlocs roots do not have a strand really
ok($location->strand, undef);
ok($location->start, 83202);
ok($location->end, 84996);

my @sublocs = $location->sub_Location();

ok(@sublocs, 2);
my $loc = shift @sublocs;
ok($loc->start, 83202);
ok($loc->end, 83329);
ok($loc->strand, -1);

$loc = shift @sublocs;
ok($loc->start, 84248);
ok($loc->end, 84996);
ok($loc->strand,1);

$seq = Bio::SeqIO->new(-format => 'GenBank' ,
                       -file=> ">".Bio::Root::IO->catfile("t","data",
                                                          "genbank.fuzzyout"));
$seq->verbose($verbosity);
ok($seq->write_seq($as));
unlink(Bio::Root::IO->catfile("t","data","genbank.fuzzyout"));

my $seqio = Bio::SeqIO->new( '-format' => 'swiss' ,
                             '-file' => Bio::Root::IO->catfile("t","data",
                                                             "swiss.dat"));

ok(defined( $seq = $seqio->next_seq));

# more tests to verify we are actually parsing correctly
ok($seq->primary_id =~ /^Bio::(Primary)?Seq/);
ok($seq->display_id, 'MA32_HUMAN');
ok($seq->length, 282);
ok($seq->division, 'HUMAN');
ok($seq->alphabet, 'protein');
my @f = $seq->all_SeqFeatures();
ok(@f, 2);
ok($f[1]->primary_tag, 'CHAIN');
ok(($f[1]->get_tag_values('description'))[0], 'COMPLEMENT COMPONENT 1, Q SUBCOMPONENT BINDING PROTEIN');

my @genenames = qw(GC1QBP HABP1 SF2P32 C1QBP);
my ($ann) = $seq->annotation->get_Annotations('gene_name');
foreach my $gn ( $ann->get_all_values() ) {
    ok ($gn, shift(@genenames));
}
ok $ann->value(-joins => [" AND "," OR "]), "GC1QBP OR HABP1 OR SF2P32 OR C1QBP";

# test for feature locations like ?..N
ok(defined( $seq = $seqio->next_seq));

ok($seq->primary_id =~ /^Bio::(Primary)?Seq/);
ok($seq->display_id, 'ACON_CAEEL');
ok($seq->length, 788);
ok($seq->division, 'CAEEL');
ok($seq->alphabet, 'protein');
ok(scalar $seq->all_SeqFeatures(), 5);

foreach my $gn ( $seq->annotation->get_Annotations('gene_name') ) {
    ok ($gn->value, 'F54H12.1');
}


# test species in swissprot -- this can be a n:n nightmare
ok ($seq = $seqio->next_seq());
my @sec_acc = $seq->get_secondary_accessions();
ok ($sec_acc[0], 'P29360');
ok ($sec_acc[1], 'Q63631');
ok ($seq->accession_number, 'P42655');
my @kw = $seq->get_keywords;
ok( $kw[0], 'Brain');
ok( $kw[1], 'Neurone');
ok ($kw[3], 'Multigene family');
ok ($seq->display_id, '143E_HUMAN');
ok ($seq->species->binomial, "Homo sapiens");
ok ($seq->species->common_name, "Human");
ok ($seq->species->ncbi_taxid, 9606);

ok ($seq = $seqio->next_seq());
ok ($seq->species->binomial, "Bos taurus");
ok ($seq->species->common_name, "Bovine");
ok ($seq->species->ncbi_taxid, 9913);

# multiple genes in swissprot
ok ($seq = $seqio->next_seq());

($ann) = $seq->annotation->get_Annotations("gene_name");
@genenames = qw(CALM1 CAM1 CALM CAM CALM2 CAM2 CAMB CALM3 CAM3 CAMC);
my $flatnames = "(CALM1 OR CAM1 OR CALM OR CAM) AND (CALM2 OR CAM2 OR CAMB) AND (CALM3 OR CAM3 OR CAMC)";

my @names = @genenames; # copy array
my @ann_names = $ann->get_all_values();

ok (scalar(@ann_names), scalar(@names));
foreach my $gn (@ann_names) {
    ok ($gn, shift(@names));
}
ok ($ann->value(-joins => [" AND "," OR "]), $flatnames);

# same entry as before, but with the new gene names format
$seqio = Bio::SeqIO->new('-format' => 'swiss' ,
                         -file => Bio::Root::IO->catfile("t","data","calm.swiss"));
ok (defined( $seq = $seqio->next_seq));
($ann) = $seq->annotation->get_Annotations("gene_name");
my @ann_names2 = $ann->get_all_values();
@names = @genenames; # copy array
ok (scalar(@ann_names2), scalar(@names));
foreach my $gn (@ann_names2) {
    ok ($gn, shift(@names));
}
ok ($ann->value(-joins => [" AND "," OR "]), $flatnames);

# test proper parsing of references
my @litrefs = $seq->annotation->get_Annotations('reference');
ok (scalar(@litrefs), 17);

my @titles = (
    '"Complete amino acid sequence of human brain calmodulin."',
    '"Multiple divergent mRNAs code for a single human calmodulin."',
    '"Molecular analysis of human and rat calmodulin complementary DNA clones. Evidence for additional active genes in these species."',
    '"Isolation and nucleotide sequence of a cDNA encoding human calmodulin."',
    '"Structure of the human CALM1 calmodulin gene and identification of two CALM1-related pseudogenes CALM1P1 and CALM1P2."',
    undef,
    '"Characterization of the human CALM2 calmodulin gene and comparison of the transcriptional activity of CALM1, CALM2 and CALM3."',
    '"Cloning of human full-length CDSs in BD Creator(TM) system donor vector."',
    '"The DNA sequence and analysis of human chromosome 14."',
    '"Generation and initial analysis of more than 15,000 full-length human and mouse cDNA sequences."',
    '"Alpha-helix nucleation by a calcium-binding peptide loop."',
    '"Solution structure of Ca(2+)-calmodulin reveals flexible hand-like properties of its domains."',
    '"Calmodulin structure refined at 1.7 A resolution."',
    '"Drug binding by calmodulin: crystal structure of a calmodulin-trifluoperazine complex."',
    '"Structural basis for the activation of anthrax adenylyl cyclase exotoxin by calmodulin."',
    '"Physiological calcium concentrations regulate calmodulin binding and catalysis of adenylyl cyclase exotoxins."',
    '"Crystal structure of a MARCKS peptide containing the calmodulin-binding domain in complex with Ca2+-calmodulin."',
);

my @locs = (
    "Biochemistry 21:2565-2569(1982).",
    "J. Biol. Chem. 263:17055-17062(1988).",
    "J. Biol. Chem. 262:16663-16670(1987).",
    "Biochem. Int. 9:177-185(1984).",
    "Eur. J. Biochem. 225:71-82(1994).",
    "Submitted (FEB-1995) to the EMBL/GenBank/DDBJ databases.",
    "Cell Calcium 23:323-338(1998).",
    "Submitted (MAY-2003) to the EMBL/GenBank/DDBJ databases.",
    "Nature 421:601-607(2003).",
    "Proc. Natl. Acad. Sci. U.S.A. 99:16899-16903(2002).",
    "Proc. Natl. Acad. Sci. U.S.A. 96:903-908(1999).",
    "Nat. Struct. Biol. 8:990-997(2001).",
    "J. Mol. Biol. 228:1177-1192(1992).",
    "Biochemistry 33:15259-15265(1994).",
    "Nature 415:396-402(2002).",
    "EMBO J. 21:6721-6732(2002).",
    "Nat. Struct. Biol. 10:226-231(2003).",
);

my @positions = (
	 undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    undef, undef,
    94, 103,
    1, 76,
    undef, undef,
    undef, undef,
    5, 148,
    1, 148,
    undef, undef,
);

foreach my $litref (@litrefs) {
    ok ($litref->title, shift(@titles));
    ok ($litref->location, shift(@locs));
    ok ($litref->start, shift(@positions));
    ok ($litref->end, shift(@positions));
}

# test dos Linefeeds in gcg parser
$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","test_badlf.gcg"),
       '-format' => 'GCG');

ok($str);
ok ( $seq = $str->next_seq());
ok(length($seq->seq) > 0 );
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $DEBUG);


$str  = new Bio::SeqIO(-format => 'genbank',
       -file   => Bio::Root::IO->catfile("t","data","AF165282.gb"),
       -verbose => $verbosity);

$seq = $str->next_seq;
@features = $seq->all_SeqFeatures();
ok(@features, 5);
ok($features[0]->start, 1);
ok($features[0]->end, 226);
$location = $features[1]->location;
ok($location->isa('Bio::Location::SplitLocationI'));
@sublocs = $location->sub_Location();
ok(@sublocs, 29);

# version and primary ID - believe it or not, this wasn't working
ok ($seq->version, 1);
ok ($seq->seq_version, 1);
ok ($seq->primary_id, "5734104");


# test embl writing of a PrimarySeq

my $primaryseq = new Bio::PrimarySeq( -seq => 'AGAGAGAGATA',
      -id  => 'myid',
      -desc => 'mydescr',
      -alphabet => 'DNA',
      -accession_number => 'myaccession');

my $embl = new Bio::SeqIO(-format => 'embl',
  -verbose => $verbosity,
  -file => ">primaryseq.embl");

ok($embl->write_seq($primaryseq));
my $scalar = "test";
eval {
    $embl->write_seq($scalar);
};
ok ($@);

unlink("primaryseq.embl");


# revcomp split location
my $gb = new Bio::SeqIO(-format => 'genbank',
                        -file   => Bio::Root::IO->catfile(qw(t data revcomp_mrna.gb)));

$seq = $gb->next_seq();

$gb = new Bio::SeqIO(-format => 'genbank',
                     -file   => ">tmp_revcomp_mrna.gb");

$gb->write_seq($seq);
undef $gb;
ok(! -z "tmp_revcomp_mrna.gb");

# INSERT DIFFING CODE HERE

unlink("tmp_revcomp_mrna.gb");


# test secondary accessions in EMBL (bug #1332)

$seqio = new Bio::SeqIO(-format =>'embl', 
                        -file => Bio::Root::IO->catfile( qw(t data ECAPAH02.embl)));
$seq = $seqio->next_seq;

ok($seq->accession_number, 'D10483');
ok($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
ok($accs[0], 'J01597');
ok($accs[-1], 'X56742');

$seqio = new Bio::SeqIO(-format => 'genbank',
                        -file=> Bio::Root::IO->catfile(qw(t data D10483.gbk)));

$seq = $seqio->next_seq;
@kw =  $seq->get_keywords;
ok(scalar @kw, 118);
ok($kw[-1], 'yabO');
@sec_acc = $seq->get_secondary_accessions();
ok(scalar @sec_acc,14);
ok($sec_acc[-1], 'X56742');


### TPA TESTS - Thanks to Richard Adams ###

### test Third Party Annotation entries in EMBL/Gb format 
#   to ensure compatability with parsers. 
###embl####

$str = new Bio::SeqIO(-format =>'embl', -file => Bio::Root::IO->catfile
      ( qw(t data BN000066-tpa.embl)));
$seq = $str->next_seq;
ok(defined $seq);
ok($seq->accession_number, 'BN000066');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'AGA000066');
ok($seq->length, 5195);
ok($seq->division, 'INV');
ok($seq->get_dates, 2);
ok($seq->keywords, 'acetylcholinesterase; achE1 gene; Third Party Annotation; TPA.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 15);

my $spec_obj = $seq->species;
ok ($spec_obj->common_name, 'African malaria mosquito');
ok ($spec_obj->species, 'gambiae');
ok ($spec_obj->genus, 'Anopheles');
ok ($spec_obj->binomial, 'Anopheles gambiae');

$ac = $seq->annotation;
my $reference =  ($ac->get_Annotations('reference') )[1];
ok ($reference->title,'"A novel acetylcholinesterase gene in mosquitoes codes for the insecticide target and is non-homologous to the ace gene in Drosophila"');
ok ($reference->authors,'Weill M., Fort P., Berthomi eu A., Dubois M.P., Pasteur N., Raymond M.');
my $cmmnt =  ($ac->get_Annotations('comment') )[0];
ok($cmmnt->text, 'see also AJ488492 for achE-1 from Kisumu strain Third Party Annotation Database: This TPA record uses Anopheles gambiae trace archive data (http://trace.ensembl.org) ');

## now genbank ##

$str = new Bio::SeqIO(-format =>'genbank', -file => Bio::Root::IO->catfile
( qw(t data BK000016-tpa.gbk)));
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->accession_number, 'BK000016');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'BK000016');
ok($seq->length, 1162);
ok($seq->division, 'ROD');
ok($seq->get_dates, 1);
ok($seq->keywords, 'Third Party Annotation; TPA');
ok($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 2);
$spec_obj = $seq->species;
ok ($spec_obj->common_name, 'Mus musculus (house mouse)');
ok ($spec_obj->species, 'musculus');
ok ($spec_obj->genus, 'Mus');
ok ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
$reference =  ($ac->get_Annotations('reference') )[0];
ok ($reference->pubmed, '11479594');
ok ($reference->medline, '21372465');

# validate that what is written is what is read

my $testfile = "testtpa.gbk";
$out = new Bio::SeqIO(-file => ">$testfile",
      -format => 'genbank');
$out->write_seq($seq);
$out->close();

$str = new Bio::SeqIO(-format =>'genbank', 
      -file => $testfile);
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
ok($seq->accession_number, 'BK000016');
ok($seq->alphabet, 'dna');
ok($seq->display_id, 'BK000016');
ok($seq->length, 1162);
ok($seq->division, 'ROD');
ok($seq->get_dates, 1);
ok($seq->keywords, 'Third Party Annotation; TPA');
ok($seq->desc, 'TPA: Mus musculus pantothenate kinase 4 mRNA, partial cds.');
ok($seq->seq_version, 1);
ok($seq->feature_count, 2);
$spec_obj = $seq->species;
ok ($spec_obj->common_name, 'Mus musculus (house mouse)');
ok ($spec_obj->species, 'musculus');
ok ($spec_obj->genus, 'Mus');
ok ($spec_obj->binomial, 'Mus musculus');
$ac = $seq->annotation;
$reference =  ($ac->get_Annotations('reference') )[0];
ok ($reference->pubmed, '11479594');
ok ($reference->medline, '21372465');

unlink($testfile);

# bug #1487

$str = new Bio::SeqIO(-verbose => $verbosity,
      -file    => Bio::Root::IO->catfile
      (qw(t data D12555.gbk)));
eval {
    $seq = $str->next_seq;
};

ok(! $@ );


# test swissprot multiple RP lines

$str = Bio::SeqIO->new(-verbose => $verbosity,
    -file    => Bio::Root::IO->catfile(qw(t data P33897) ));

ok($seq = $str->next_seq);
my @refs = $seq->annotation->get_Annotations('reference');
ok( @refs, 23);
ok($refs[20]->rp, 'VARIANTS X-ALD LEU-98; ASP-99; GLU-217; GLN-518; ASP-608; ILE-633 AND PRO-660, AND VARIANT THR-13.');


# test bug #1673 , RDB-II genbank files
$str = Bio::SeqIO->new(-format => 'genbank', 
                       -file => Bio::Root::IO->catfile(qw(t data Mcjanrna_rdbII.gbk) )
		      );
ok($seq = $str->next_seq);
@refs = $seq->annotation->get_Annotations('reference');
ok(@refs, 1);
ok($seq->display_id,'Mc.janrrnA');
ok($seq->molecule ,'RNA');


# test BSML-SAX
unless( $NOSAX ){
	$str = Bio::SeqIO->new(-format => 'bsml_sax', 
			       -file => Bio::Root::IO->catfile
			       (qw(t data U83300.bsml) )
			       );
	ok($seq = $str->next_seq);
	@refs = $seq->annotation->get_Annotations('reference');
	ok(@refs, 2);
	ok($seq->display_id,'MIVN83300');
	ok($seq->molecule ,'dna');
	ok(! $seq->is_circular);
	ok($seq->get_dates,2);
	ok($seq->accession_number, 'U83300');
	ok($seq->seq_version,1);
	my @feats = $seq->get_SeqFeatures;
	ok(@feats, 2);
	ok($feats[1]->start, 1);
	ok($feats[1]->end, 946);
	ok($feats[1]->get_tag_values('db_xref'), 3);
	ok($seq->annotation->get_Annotations('reference'),2);
	ok($seq->annotation->get_Annotations('dblink'),2);
}


# bug 1571
{
    my $ent = Bio::SeqIO->new(-format => 'embl',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data test.embl2sq)));
    ok($ent->next_seq->length,4877);
}
