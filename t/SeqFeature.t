# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use constant NUMTESTS => 233;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;


my $skipdbtests;
my $skip_all;
BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    
    plan tests => NUMTESTS;
    eval { 
	require IO::String; 
	require LWP::UserAgent;
	require HTTP::Request::Common;
	require Bio::DB::GenBank;
    };
    if( $@ ) {
	print STDERR "IO::String, LWP::UserAgent or HTTP::Request not installed - skipping DB tests...\n";
	$skipdbtests = 1;
    } else {
	$skipdbtests = 0;
    }
    eval {
	require URI::Escape;
    };
    if( $@ ) {
	print STDERR "URI::Escape not installed, so Bio::SeqFeature::Annotated not usable - skipping all tests...\n";
	$skip_all = 1;
    }
}

use_ok('Bio::Seq');
use_ok('Bio::SeqIO');
use_ok('Bio::SeqFeature::Generic');
use_ok('Bio::SeqFeature::FeaturePair');
use_ok('Bio::SeqFeature::SimilarityPair');
use_ok('Bio::SeqFeature::Computation');
use_ok('Bio::SeqFeature::Gene::Transcript');
use_ok('Bio::SeqFeature::Gene::UTR');
use_ok('Bio::SeqFeature::Gene::Exon');
use_ok('Bio::SeqFeature::Gene::Poly_A_site');
use_ok('Bio::SeqFeature::Gene::GeneStructure');
use_ok('Bio::Location::Fuzzy');

# predeclare variables for strict
my ($feat,$str,$feat2,$pair,$comp_obj1,$comp_obj2,@sft); 


my $verbose = 0;

$feat = new Bio::SeqFeature::Generic ( -start => 40,
				       -end => 80,
				       -strand => 1,
				       -primary => 'exon',
				       -source  => 'internal',
				       -tag => {
					   silly => 20,
					   new => 1
					   }
				       );

is $feat->start, 40, 'start of feature location';
is $feat->end, 80, 'end of feature location';
is $feat->primary_tag, 'exon', 'primary tag';
is $feat->source_tag, 'internal', 'source tag';

$str = $feat->gff_string() || ""; # placate -w

$pair = new Bio::SeqFeature::FeaturePair();
ok defined $pair;

$feat2 = new Bio::SeqFeature::Generic ( -start => 400,
				       -end => 440,
				       -strand => 1,
				       -primary => 'other',
				       -source  => 'program_a',
				       -tag => {
					   silly => 20,
					   new => 1
					   }
				       );

ok defined $feat2;
$pair->feature1($feat);
$pair->feature2($feat2);

is $pair->feature1, $feat, 'feature1 of pair stored';
is $pair->feature2, $feat2, 'feature2 of pair stored';
is $pair->start, 40, 'feature start';
is $pair->end, 80, 'feature end';
is $pair->primary_tag, 'exon', 'primary tag';
is $pair->source_tag, 'internal', 'source tag';
is $pair->hstart, 400, 'hstart';
is $pair->hend, 440, 'hend';
is $pair->hprimary_tag, 'other', 'hprimary tag';
is $pair->hsource_tag, 'program_a', 'hsource tag';

$pair->invert;
is $pair->end, 440, 'inverted end';

# Test attaching a SeqFeature::Generic to a Bio::Seq
{
    # Make the parent sequence object
    my $seq = Bio::Seq->new(
			    -seq          => 'aaaaggggtttt',
			    -display_id   => 'test',
			    -alphabet     => 'dna',
			    );
    
    # Make a SeqFeature
    my $sf1 = Bio::SeqFeature::Generic->new(
					    -start    => 4,
					    -end      => 9,
					    -strand   => 1,
					    );
    
    # Add the SeqFeature to the parent
    ok ($seq->add_SeqFeature($sf1));
    
    # Test that it gives the correct sequence
    my $sf_seq1 = $sf1->seq->seq;
    is $sf_seq1, 'aggggt', 'seq string';
    is $sf1->end,9, 'sf1 end';
    is $sf1->start,4, 'sf1 start';

    # Make a second seqfeature on the opposite strand
    my $sf2 = Bio::SeqFeature::Generic->new(
					    -start    => 4,
					    -end      => 9,
					    -strand   => -1,
					    );
    
    # This time add the PrimarySeq to the seqfeature
    # before adding it to the parent
    ok ($sf2->attach_seq($seq->primary_seq));
    $seq->add_SeqFeature($sf2);
    
    # Test again that we have the correct sequence
    my $sf_seq2 = $sf2->seq->seq;
    is $sf_seq2, 'acccct', 'sf2';
}

#Do some tests for computation.pm

ok defined ( $comp_obj1 = Bio::SeqFeature::Computation->new('-start' => 1,
							    '-end'   => 10) );
is($comp_obj1->computation_id(332),332, 'computation id');
ok( $comp_obj1->add_score_value('P', 33), 'score value');
{
    $comp_obj2 = Bio::SeqFeature::Computation->new('-start' => 2,
						   '-end'   => 10);
    ok ($comp_obj1->add_sub_SeqFeature($comp_obj2, 'exon') );
    ok (@sft = $comp_obj1->all_sub_SeqFeature_types() );
    is($sft[0], 'exon', 'sft[0] is exon');
}

ok defined ( $comp_obj1 = new Bio::SeqFeature::Computation 
	     (
	      -start => 10, -end => 100,
	      -strand => -1, -primary => 'repeat',
	      -program_name => 'GeneMark',
	      -program_date => '12-5-2000',
	      -program_version => 'x.y',
	      -database_name => 'Arabidopsis',
	      -database_date => '12-dec-2000',
	      -computation_id => 2231,
	      -score    => { no_score => 334 } )
	     );

is ( $comp_obj1->computation_id, 2231, 'computation id' );
ok ( $comp_obj1->add_score_value('P', 33) );
is ( ($comp_obj1->each_score_value('no_score'))[0], '334', 'score value');

# some tests for bug #947

my $sfeat = new Bio::SeqFeature::Generic(-primary => 'test');

$sfeat->add_sub_SeqFeature(new Bio::SeqFeature::Generic(-start => 2,
							-end   => 4,
							-primary => 'sub1'),
			   'EXPAND');

$sfeat->add_sub_SeqFeature(new Bio::SeqFeature::Generic(-start => 6,
							-end   => 8,
							-primary => 'sub2'),
			   'EXPAND');

is $sfeat->start, 2, 'sfeat start for EXPAND-ED feature (bug #947)';
is $sfeat->end, 8, 'sfeat end for EXPAND-ED feature (bug #947)';

# some tests to see if we can set a feature to start at 0
$sfeat = new Bio::SeqFeature::Generic(-start => 0, -end => 0 );

ok(defined $sfeat->start);
is($sfeat->start,0, 'can create feature starting and ending at 0');
ok(defined $sfeat->end);
is($sfeat->end,0,'can create feature starting and ending at 0');


# tests for Bio::SeqFeature::Gene::* objects
# using information from acc: AB077698 as a guide

ok my $seqio = new Bio::SeqIO(-format => 'genbank',
			 -file   => Bio::Root::IO->catfile("t","data","AB077698.gb"));
ok my $geneseq = $seqio->next_seq();

ok my $gene = new Bio::SeqFeature::Gene::GeneStructure(-primary => 'gene',
						    -start   => 1,
						       -end     => 2701,
						       -strand  => 1);

ok my $transcript = new Bio::SeqFeature::Gene::Transcript(-primary => 'CDS',
							  -start   => 80,
							  -end     => 1144,
							  -tag     => { 
							      'gene' => "CHCR",
							      'note' => "Cys3His CCG1-Required Encoded on BAC clone RP5-842K24 (AL050310) The human CHCR (Cys3His CCG1-Required) protein is highly related to EXP/MBNL (Y13829, NM_021038, AF401998) and MBLL (NM_005757,AF061261), which together comprise the human Muscleblind family",
							      'codon_start' => 1,
							      'protein_id'  => 'BAB85648.1',
							  });

ok my $poly_A_site1 = new Bio::SeqFeature::Gene::Poly_A_site
    (-primary => 'polyA_site',
     -start => 2660,
     -end   => 2660,
     -tag   => { 
	 'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#2 used by CHCR EST clone DKFZp434G2222 (AL133625)"
	 });

ok my $poly_A_site2 = new Bio::SeqFeature::Gene::Poly_A_site
    (-primary => 'polyA_site',
     -start => 1606,
     -end   => 1606,
     -tag   => { 
	 'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#1 used by CHCR EST clone PLACE1010202 (AK002178)",
     });

ok my $fiveprimeUTR = new Bio::SeqFeature::Gene::UTR(-primary => "utr5prime");
ok $fiveprimeUTR->location(new Bio::Location::Fuzzy(-start => "<1",
						    -end   => 79));
ok my $threeprimeUTR = new Bio::SeqFeature::Gene::UTR(-primary => "utr3prime",
						      -start   => 1145,
						      -end     => 2659);

# Did a quick est2genome against genomic DNA (this is on Chr X) to
# get the gene structure by hand since it is not in the file
# --Jason

ok my $exon1 = new Bio::SeqFeature::Gene::Exon(-primary => 'exon',
					       -start => 80,
					       -end   => 177);
ok $geneseq->add_SeqFeature($exon1);

ok $geneseq->add_SeqFeature($fiveprimeUTR);
ok $geneseq->add_SeqFeature($threeprimeUTR);
ok $geneseq->add_SeqFeature($poly_A_site1);
ok $geneseq->add_SeqFeature($poly_A_site2);

ok $transcript->add_utr($fiveprimeUTR, 'utr5prime');
ok $transcript->add_utr($threeprimeUTR, 'utr3prime');

ok $transcript->add_exon($exon1);

# API only supports a single poly-A site per transcript at this point 
$transcript->poly_A_site($poly_A_site2);
$geneseq->add_SeqFeature($transcript);
$gene->add_transcript($transcript);
$geneseq->add_SeqFeature($gene);

my ($t) = $gene->transcripts(); # get 1st transcript
ok(defined $t); 
is($t->mrna->length, 1693, 'mRNA spliced length');
is($gene->utrs, 2, 'has 2 UTRs');



# Test for bug when Locations are not created explicitly

my $feat1 = new Bio::SeqFeature::Generic(-start => 1,
					 -end   => 15,
					 -strand=> 1);

$feat2 = new Bio::SeqFeature::Generic(-start => 10,
					 -end   => 25,
					 -strand=> 1);

my $overlap = $feat1->location->union($feat2->location);
is($overlap->start, 1);
is($overlap->end,   25);

my $intersect = $feat1->location->intersection($feat2->location);
is($intersect->start, 10);
is($intersect->end,   15);


# now let's test spliced_seq

isa_ok(  $seqio = new Bio::SeqIO(-file => Bio::Root::IO->catfile
				(qw(t data AY095303S1.gbk)),
				 -format  => 'genbank'), "Bio::SeqIO");

isa_ok( $geneseq = $seqio->next_seq(), 'Bio::Seq');
my ($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
my $db;

if( $skipdbtests ) {
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',4);
} else {
    $db = new Bio::DB::GenBank(-verbose=> $DEBUG);
    $CDS->verbose(-1);
    my $cdsseq = $CDS->spliced_seq(-db => $db,-nosort => 1);
    
    is($cdsseq->subseq(1,76),
       'ATGCAGCCATACGCTTCCGTGAGCGGGCGATGTCTATCTAGACCAGATGCATTGCATGTGATACCGTTTGGGCGAC');
    is($cdsseq->translate->subseq(1,100), 
       'MQPYASVSGRCLSRPDALHVIPFGRPLQAIAGRRFVRCFAKGGQPGDKKKLNVTDKLRLGNTPPTLDVLKAPRPTDAPSAIDDAPSTSGLGLGGGVASPR');
    # test what happens without 
    $cdsseq = $CDS->spliced_seq(-db => $db,-nosort => 1);    
    is($cdsseq->subseq(1,76), 
       'ATGCAGCCATACGCTTCCGTGAGCGGGCGATGTCTATCTAGACCAGATGCATTGCATGTGATACCGTTTGGGCGAC');
    is($cdsseq->translate->subseq(1,100), 
       'MQPYASVSGRCLSRPDALHVIPFGRPLQAIAGRRFVRCFAKGGQPGDKKKLNVTDKLRLGNTPPTLDVLKAPRPTDAPSAIDDAPSTSGLGLGGGVASPR');
    
} 

isa_ok(  $seqio = new Bio::SeqIO(-file => Bio::Root::IO->catfile
				(qw(t data AF032047.gbk)),
				-format  => 'genbank'), 'Bio::SeqIO');
isa_ok $geneseq = $seqio->next_seq(), 'Bio::Seq';
($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
if ($skipdbtests ) { 
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',2);
} else {
    my $cdsseq = $CDS->spliced_seq( -db => $db, -nosort => 1);
    is($cdsseq->subseq(1,70), 'ATGGCTCGCTTCGTGGTGGTAGCCCTGCTCGCGCTACTCTCTCTGTCTGGCCTGGAGGCTATCCAGCATG');
    is($cdsseq->translate->seq, 'MARFVVVALLALLSLSGLEAIQHAPKIQVYSRHPAENGKPNFLNCYVSGFHPSDIEVDLLKNGKKIEKVEHSDLSFSKDWSFYLLYYTEFTPNEKDEYACRVSHVTFPTPKTVKWDRTM*');
}


# trans-spliced 

isa_ok( $seqio = Bio::SeqIO->new(-format => 'genbank',
				 -file   => 
				 Bio::Root::IO->catfile(qw(t data NC_001284.gbk))), 
	'Bio::SeqIO');
my $genome = $seqio->next_seq;

foreach my $cds (grep { $_->primary_tag eq 'CDS' } $genome->get_SeqFeatures) {
   my $spliced = $cds->spliced_seq(-nosort => 1)->translate->seq;
   chop($spliced); # remove stop codon
   is($spliced,($cds->get_tag_values('translation'))[0],
      'spliced seq translation matches expected');
}

# spliced_seq phase 
my $seqin = Bio::SeqIO->new(-format => 'fasta',
                            -file   =>
                            Bio::Root::IO->catfile(qw(t data sbay_c127.fas)));

my $seq = $seqin->next_seq;

my $sf = Bio::SeqFeature::Generic->new(-verbose => -1,
                                       -start => 263,
                                       -end => 721,
                                       -strand => 1,
                                       -primary => 'splicedgene');

$sf->attach_seq($seq);

my %phase_check = (
    'TTCAATGACT' => 'FNDFYSMGKS',
    'TCAATGACTT' => 'SMTSIPWVNQ',
    'GTTCAATGAC' => 'VQ*LLFHG*I',
);

for my $phase (-1..3) {
    my $sfseq = $sf->spliced_seq(-phase => $phase);
    ok exists $phase_check{$sfseq->subseq(1,10)};
    is ($sfseq->translate->subseq(1,10), $phase_check{$sfseq->subseq(1,10)}, 'phase check');
}


SKIP: {
# use_ok('Bio::SeqFeature::Annotated');
    eval { require Bio::SeqFeature::Annotated};
    
    my $sfa = Bio::SeqFeature::Annotated->new(-start => 1,
					      -end => 5,
					      -strand => "+",
					      -frame => 2,
					      -phase => 2,
					      -score => 12,
					      -display_name => 'test.annot',
					      -seq_id => 'test.displayname' );
    
    ok (defined $sfa);
    my $loc = $sfa->location;
    ok $loc->isa("Bio::Location::Simple");    
    ok $sfa->display_name eq 'test.annot';

#test bsfa::from_feature
    my $sfg = Bio::SeqFeature::Generic->new ( -start => 400,
					      -end => 440,
					      -strand => 1,
					      -primary => 'nucleotide_motif',
					      -source  => 'program_a',
					      -tag => {
						  silly => 20,
						  new => 1
						  }
					      );
    my $sfa2;
    $sfa2 = Bio::SeqFeature::Annotated->new(-feature => $sfg);
    
    is $sfa2->type->name,'nucleotide_motif';
    is $sfa2->primary_tag, 'nucleotide_motif';
    is $sfa2->source,'program_a';
    is $sfa2->strand,1;
    is $sfa2->start,400;
    is $sfa2->end,440;
    is $sfa2->get_Annotations('silly')->value,20;
    is $sfa2->get_Annotations('new')->value,1;
    
    my $sfa3 = Bio::SeqFeature::Annotated->new( -start => 1,
						-end => 5,
						-strand => "+",
						-frame => 2,
						-phase => 2,
						-score => 12,
						-display_name => 'test.annot',
						-seq_id => 'test.displayname' );
    $sfa3->from_feature($sfg);
    
    
    is $sfa3->type->name,'nucleotide_motif', 'type->name';
    is $sfa3->primary_tag, 'nucleotide_motif', 'primary_tag';
    is $sfa3->source,'program_a';
    is $sfa3->strand,1;
    is $sfa3->start,400;
    is $sfa3->end,440;
    is $sfa3->get_Annotations('silly')->value,20;
    is $sfa3->get_Annotations('new')->value,1;
}
