# -*-Perl-*- mode (to keep my emacs happy)

use strict;
use vars qw($DEBUG);
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 90;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;
use Bio::Root::IO;
use Bio::Annotation;

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
ok($as->primary_id, 'ROA1');
ok($as->length, 371);
ok($as->alphabet, 'protein');
ok($as->division, 'HUMAN');
ok(scalar $as->all_SeqFeatures(), 16);

($ent, $seq, $out) = undef;

$ent = Bio::SeqIO->new( '-file' => Bio::Root::IO->catfile("t","data","test.embl"),
			'-format' => 'embl');

$seq = $ent->next_seq();

ok(defined $seq->seq(), 1, 
   'failure to read Embl with ^ location and badly split double quotes');

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
while($seq = $stream->next_seq()) {
    $seqnum++;
    if($seqnum == 3) {
	ok $seq->display_id(), "HUMBDNF";
	$species = $seq->species();
	@cl = $species->classification();
	ok( $species->binomial(), "Homo sapiens", 
	    'species parsing incorrect for genbank');
	ok( $cl[3] ne $species->genus(), 1, 
	    'genus duplicated in genbank parsing');
    }
    $lasts = $seq;
}
ok $lasts->display_id(), "HUMBETGLOA";
$stream->close();
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

ok($lastfeature->strand, -1);
my $location = $lastfeature->location;
ok($location->strand, -1);
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
ok($seq->primary_id, 'MA32');
ok($seq->display_id, 'MA32_HUMAN');
ok($seq->length, 282);
ok($seq->division, 'HUMAN');
ok($seq->alphabet, 'protein');
ok(scalar $seq->all_SeqFeatures(), 2);

my $seen = 0;
foreach my $gn ( $seq->annotation->get_Annotations('gene_name') ) {
    if( $gn->value =~ /SF2/ ) {
	$seen = 1;
    }	       
}

ok $seen;
# test for feature locations like ?..N
ok(defined( $seq = $seqio->next_seq));

ok($seq->primary_id, 'ACON');
ok($seq->display_id, 'ACON_CAEEL');
ok($seq->length, 788);
ok($seq->division, 'CAEEL');
ok($seq->alphabet, 'protein');
ok(scalar $seq->all_SeqFeatures(), 5);

$seen = 0;
foreach my $gn ( $seq->annotation->get_Annotations('gene_name') ) {
    if( $gn->value =~ /F54H12/ ) {
	$seen = 1;
    }	       
}
ok($seen);

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
				      -alphabet => 'DNA',
				      -accession_number => 'myaccession');
my $embl = new Bio::SeqIO(-format => 'embl' );

#ok($embl->write_seq($primaryseq));
