# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
use vars qw($DEBUG $TESTCOUNT);
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 161;
    plan tests => $TESTCOUNT;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use Bio::Root::IO;
use Bio::Annotation::Collection;

ok(1);

my $verbosity = -1;   # Set to -1 for release version, so warnings aren't printed

my ($str, $seq,$ast,$temp,$mf,$ent,$out); # predeclare variables for strict
$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","test.fasta"), 
		       '-format' => 'Fasta');
ok $str;

ok (defined($seq = $str->next_seq()));

print "Sequence 1 of 2 from fasta stream:\n", $seq->seq, "\n" if ( $DEBUG);

ok($seq->id, 'roa1_drome');
ok $seq->length, 358;


$str = Bio::SeqIO->new(-file=> Bio::Root::IO->catfile("t","data","test.raw"), '-format' => 'Raw');

ok $str;

ok ($seq = $str->next_seq());
print "Sequence 1 of 2 from Raw stream:\n", $seq->seq, "\n\n" if( $DEBUG);

ok ($seq = $str->next_seq());
    
print "Sequence 2 of 2 from Raw stream:\n", $seq->seq, $seq->seq, "\n" 
    if( $DEBUG);



$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","test.gcg"), 
		       '-format' => 'GCG');

ok $str;

ok ( $seq = $str->next_seq());
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $DEBUG);


$str = Bio::SeqIO->new('-file'=> ">".Bio::Root::IO->catfile("t","data","gcg.out"), 
		       '-format' => 'GCG');

$str->write_seq($seq);
ok(1);
unlink(Bio::Root::IO->catfile("t","data","gcg.out"));


$str = Bio::SeqIO->new( '-file'=> Bio::Root::IO->catfile("t","data","test.genbank"), 
			'-format' => 'GenBank');

ok $str;
$str->verbose($verbosity);

ok ( $seq = $str->next_seq() );
print "Sequence 1 of 1 from GenBank stream:\n", $seq->seq, "\n" if( $DEBUG);


my $strout = Bio::SeqIO->new('-file'=> ">".Bio::Root::IO->catfile("t","data","genbank.out"), 
			     '-format' => 'GenBank');
while( $seq ) {
    $strout->write_seq($seq);
    $seq = $str->next_seq();
}
undef $strout;
unlink(Bio::Root::IO->catfile("t","data","genbank.out"));

ok(1);

$str = undef;


$ast = Bio::SeqIO->new( '-format' => 'embl' , 
			'-file' => Bio::Root::IO->catfile("t","data","roa1.dat"));
$ast->verbose($verbosity);
my $as = $ast->next_seq();
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

$ast = Bio::SeqIO->new( '-format' => 'GenBank' , 
			'-file' => Bio::Root::IO->catfile("t","data","roa1.genbank"));
$ast->verbose($verbosity);
$as = $ast->next_seq();
ok $as->molecule, 'mRNA';
ok $as->alphabet, 'dna';

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
$ast = Bio::SeqIO->new( '-verbosity' => $verbosity,
			'-format' => 'swiss' , 
			'-file' => Bio::Root::IO->catfile("t","data","roa1.swiss"));
$as = $ast->next_seq();

ok defined $as->seq;
ok($as->id, 'ROA1_HUMAN', "id is ".$as->id);
#ok($as->primary_id, 'ROA1');
skip($as->primary_id =~ /^Bio::Seq::/, $as->primary_id, 'ROA1');
ok($as->length, 371);
ok($as->alphabet, 'protein');
ok($as->division, 'HUMAN');
ok(scalar $as->all_SeqFeatures(), 16);

ok(scalar $as->annotation->get_Annotations('reference'), 11);

($ent, $seq, $out,$as) = undef;

$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data","test.embl"),
			'-format' => 'embl');

$seq = $ent->next_seq();

ok(defined $seq->seq(), 1, 
   'failure to read Embl with ^ location and badly split double quotes');
ok(scalar $seq->annotation->get_Annotations('reference'), 3);
$out = Bio::SeqIO->new('-file'=> ">". Bio::Root::IO->catfile("t","data","embl.out"), 
		       '-format' => 'embl');

ok($out->write_seq($seq),1,
   'failure to write Embl format with ^ < and > locations');

unlink(Bio::Root::IO->catfile("t","data","embl.out"));

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

my $stream = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","test.genbank"),
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

$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data","test.embl"), 
			'-format' => 'embl');
$ent->verbose($verbosity);
$seq = $ent->next_seq();
$species = $seq->species();
@cl = $species->classification();
ok( $cl[3] ne $species->genus(), 1, 'genus duplicated in EMBL parsing');
$ent->close();

$seq = Bio::SeqIO->new( '-format' => 'GenBank' , 
			-file => Bio::Root::IO->catfile("t","data","testfuzzy.genbank"));
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

$seq = Bio::SeqIO->new( '-format' => 'GenBank' , 
			-file => ">".Bio::Root::IO->catfile("t","data","genbank.fuzzyout"));
$seq->verbose($verbosity);
ok($seq->write_seq($as));
unlink(Bio::Root::IO->catfile("t","data","genbank.fuzzyout"));

my $seqio = Bio::SeqIO->new( '-format' => 'swiss' ,
			  -file => Bio::Root::IO->catfile("t","data","swiss.dat"));

ok(defined( $seq = $seqio->next_seq));

# more tests to verify we are actually parsing correctly
skip($seq->primary_id =~ /^Bio::Seq/, $seq->primary_id, 'MA32');
ok($seq->display_id, 'MA32_HUMAN');
ok($seq->length, 282);
ok($seq->division, 'HUMAN');
ok($seq->alphabet, 'protein');
ok(scalar $seq->all_SeqFeatures(), 2);

my @genenames = qw(GC1QBP HABP1 SF2P32 C1QBP);
my ($ann) = $seq->annotation->get_Annotations('gene_name');
foreach my $gn ( $ann->get_all_values() ) {
    ok ($gn, shift(@genenames));
}
ok $ann->value(-joins => [" AND "," OR "]), "GC1QBP OR HABP1 OR SF2P32 OR C1QBP";

# test for feature locations like ?..N
ok(defined( $seq = $seqio->next_seq));

skip($seq->primary_id =~ /^Bio::Seq/, $seq->primary_id, 'ACON');
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
foreach my $gn ( $ann->get_all_values() ) {
    ok ($gn, shift(@genenames));
}
ok $ann->value(-joins => [" AND "," OR "]), "(CALM1 OR CAM1 OR CALM OR CAM) AND (CALM2 OR CAM2 OR CAMB) AND (CALM3 OR CAM3 OR CAMC)";

# test dos Linefeeds in gcg parser
$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","test_badlf.gcg"), 
		       '-format' => 'GCG');

ok($str);
ok ( $seq = $str->next_seq());
ok(length($seq->seq) > 0 );
print "Sequence 1 of 1 from GCG stream:\n", $seq->seq, "\n" if( $DEBUG);


$str  = new Bio::SeqIO(-format => 'genbank',
		       -file   => Bio::Root::IO->catfile("t","data",
							 "AF165282.gb"),
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

# PIR testing

$str = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile("t","data","seqfile.pir"), 
		       '-format' => 'pir');
ok $str;
$strout = new Bio::SeqIO(-format => 'pir', -fh => \*STDOUT);

while( $seq = $str->next_seq()) {
    ok($seq->id, qr /^[PF]1/ );
    ok($seq->length > 1);
    $strout->write_seq($seq) if( $verbosity > 0);
}

# test embl writing of a PrimarySeq

my $primaryseq = new Bio::PrimarySeq( -seq => 'AGAGAGAGATA',
				      -id  => 'myid',
				      -desc => 'mydescr',
				      -alphabet => 'DNA',
				      -accession_number => 'myaccession');

my $embl = new Bio::SeqIO(-format => 'embl', 
			  -verbose => $verbosity -1,
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

$seqio = new Bio::SeqIO(-format =>'embl', -file => Bio::Root::IO->catfile
			( qw(t data ECAPAH02.embl)));
$seq = $seqio->next_seq;

ok($seq->accession_number, 'D10483');
ok($seq->seq_version, 2);
my @accs = $seq->get_secondary_accessions();
ok($accs[0], 'J01597');
ok($accs[-1], 'X56742');

$seqio = new Bio::SeqIO(-format => 'genbank',
			-file   => Bio::Root::IO->catfile(qw(t data 
							     D10483.gbk)));

$seq = $seqio->next_seq;
my @kw =  $seq->get_keywords;
ok(scalar @kw, 118);
ok($kw[-1], 'yabO');
my @sec_acc = $seq->get_secondary_accessions();
ok(scalar @sec_acc,23);
ok($sec_acc[-1], 'X56742');
