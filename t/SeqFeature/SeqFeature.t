# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 245);
	
	use_ok('Bio::Seq');
	use_ok('Bio::SeqIO');
	use_ok('Bio::SeqFeature::Generic');
	use_ok('Bio::SeqFeature::FeaturePair');
	use_ok('Bio::SeqFeature::Computation');
	use_ok('Bio::SeqFeature::Gene::Transcript');
	use_ok('Bio::SeqFeature::Gene::UTR');
	use_ok('Bio::SeqFeature::Gene::Exon');
	use_ok('Bio::SeqFeature::Gene::Poly_A_site');
	use_ok('Bio::SeqFeature::Gene::GeneStructure');
	use_ok('Bio::Location::Fuzzy');
}

# predeclare variables for strict
my ($feat,$str,$feat2,$pair,$comp_obj1,$comp_obj2,@sft); 

my $DEBUG = test_debug();

$feat = Bio::SeqFeature::Generic->new( -start => 40,
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

$pair = Bio::SeqFeature::FeaturePair->new();
ok defined $pair;

$feat2 = Bio::SeqFeature::Generic->new( -start => 400,
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

ok defined ( $comp_obj1 = Bio::SeqFeature::Computation->new 
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

my $sfeat = Bio::SeqFeature::Generic->new(-primary => 'test');

$sfeat->add_sub_SeqFeature(Bio::SeqFeature::Generic->new(-start => 2,
							-end   => 4,
							-primary => 'sub1'),
			   'EXPAND');

$sfeat->add_sub_SeqFeature(Bio::SeqFeature::Generic->new(-start => 6,
							-end   => 8,
							-primary => 'sub2'),
			   'EXPAND');

is $sfeat->start, 2, 'sfeat start for EXPAND-ED feature (bug #947)';
is $sfeat->end, 8, 'sfeat end for EXPAND-ED feature (bug #947)';

# some tests to see if we can set a feature to start at 0
$sfeat = Bio::SeqFeature::Generic->new(-start => 0, -end => 0 );

ok(defined $sfeat->start);
is($sfeat->start,0, 'can create feature starting and ending at 0');
ok(defined $sfeat->end);
is($sfeat->end,0,'can create feature starting and ending at 0');


# tests for Bio::SeqFeature::Gene::* objects
# using information from acc: AB077698 as a guide

ok my $seqio = Bio::SeqIO->new(-format => 'genbank',
			 -file   => test_input_file('AB077698.gb'));
ok my $geneseq = $seqio->next_seq();

ok my $gene = Bio::SeqFeature::Gene::GeneStructure->new(-primary => 'gene',
						    -start   => 1,
						       -end     => 2701,
						       -strand  => 1);

ok my $transcript = Bio::SeqFeature::Gene::Transcript->new(-primary => 'CDS',
							  -start   => 80,
							  -end     => 1144,
							  -tag     => { 
							      'gene' => "CHCR",
							      'note' => "Cys3His CCG1-Required Encoded on BAC clone RP5-842K24 (AL050310) The human CHCR (Cys3His CCG1-Required) protein is highly related to EXP/MBNL (Y13829, NM_021038, AF401998) and MBLL (NM_005757,AF061261), which together comprise the human Muscleblind family",
							      'codon_start' => 1,
							      'protein_id'  => 'BAB85648.1',
							  });

ok my $poly_A_site1 = Bio::SeqFeature::Gene::Poly_A_site->new
    (-primary => 'polyA_site',
     -start => 2660,
     -end   => 2660,
     -tag   => { 
	 'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#2 used by CHCR EST clone DKFZp434G2222 (AL133625)"
	 });

ok my $poly_A_site2 = Bio::SeqFeature::Gene::Poly_A_site->new
    (-primary => 'polyA_site',
     -start => 1606,
     -end   => 1606,
     -tag   => { 
	 'note' => "Encoded on BAC clone RP5-842K24 (AL050310); PolyA_site#1 used by CHCR EST clone PLACE1010202 (AK002178)",
     });

ok my $fiveprimeUTR = Bio::SeqFeature::Gene::UTR->new(-primary => "utr5prime");
ok $fiveprimeUTR->location(Bio::Location::Fuzzy->new(-start => "<1",
						    -end   => 79));
ok my $threeprimeUTR = Bio::SeqFeature::Gene::UTR->new(-primary => "utr3prime",
						      -start   => 1145,
						      -end     => 2659);

# Did a quick est2genome against genomic DNA (this is on Chr X) to
# get the gene structure by hand since it is not in the file
# --Jason

ok my $exon1 = Bio::SeqFeature::Gene::Exon->new(-primary => 'exon',
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

my $feat1 = Bio::SeqFeature::Generic->new(-start => 1,
					 -end   => 15,
					 -strand=> 1);

$feat2 = Bio::SeqFeature::Generic->new(-start => 10,
					 -end   => 25,
					 -strand=> 1);

my $overlap = $feat1->location->union($feat2->location);
is($overlap->start, 1);
is($overlap->end,   25);

my $intersect = $feat1->location->intersection($feat2->location);
is($intersect->start, 10);
is($intersect->end,   15);


# now let's test spliced_seq

isa_ok(  $seqio = Bio::SeqIO->new(-file => test_input_file('AY095303S1.gbk'),
				 -format  => 'genbank'), "Bio::SeqIO");

isa_ok( $geneseq = $seqio->next_seq(), 'Bio::Seq');
my ($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
my $db;

SKIP: {
	test_skip(-tests => 5,
			  -requires_modules => [qw(IO::String
									   LWP::UserAgent
									   HTTP::Request::Common)],
			  -requires_networking => 1);
	
	use_ok('Bio::DB::GenBank');
    
    $db = Bio::DB::GenBank->new(-verbose=> $DEBUG);
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

isa_ok(  $seqio = Bio::SeqIO->new(-file => test_input_file('AF032047.gbk'),
				-format  => 'genbank'), 'Bio::SeqIO');
isa_ok $geneseq = $seqio->next_seq(), 'Bio::Seq';
($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
SKIP: { 
    test_skip(-tests => 2,
			  -requires_modules => [qw(IO::String
									   LWP::UserAgent
									   HTTP::Request::Common)],
			  -requires_networking => 1);
	
    my $cdsseq = $CDS->spliced_seq( -db => $db, -nosort => 1);
    is($cdsseq->subseq(1,70), 'ATGGCTCGCTTCGTGGTGGTAGCCCTGCTCGCGCTACTCTCTCTGTCTGGCCTGGAGGCTATCCAGCATG');
    is($cdsseq->translate->seq, 'MARFVVVALLALLSLSGLEAIQHAPKIQVYSRHPAENGKPNFLNCYVSGFHPSDIEVDLLKNGKKIEKVEHSDLSFSKDWSFYLLYYTEFTPNEKDEYACRVSHVTFPTPKTVKWDRTM*');
}


# trans-spliced 

isa_ok( $seqio = Bio::SeqIO->new(-format => 'genbank',
				 -file   => test_input_file('NC_001284.gbk')), 
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
                            -file   => test_input_file('sbay_c127.fas'));

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

# tags
$sf->add_tag_value('note','n1');
$sf->add_tag_value('note','n2');
$sf->add_tag_value('comment','c1');
is_deeply( [sort $sf->get_all_tags()], [sort qw(note comment)] , 'tags found');
is_deeply( [sort $sf->get_tagset_values('note')], [sort qw(n1 n2)] , 'get_tagset_values tag values found');
is_deeply( [sort $sf->get_tagset_values(qw(note comment))], [sort qw(c1 n1 n2)] , 'get_tagset_values tag values for multiple tags found');
lives_ok { 
  is_deeply( [sort $sf->get_tag_values('note')], [sort qw(n1 n2)] , 'get_tag_values tag values found');
} 'get_tag_values lives with tag';
lives_ok { 
  is_deeply( [$sf->get_tagset_values('notag') ], [], 'get_tagset_values no tag values found');
} 'get_tagset_values lives with no tag';
throws_ok { $sf->get_tag_values('notag') } qr/tag value that does not exist/, 'get_tag_values throws with no tag';

# circular sequence SeqFeature tests
$seqin = Bio::SeqIO->new(-format => 'genbank',
                         -file   => test_input_file('PX1CG.gb'));

$seq = $seqin->next_seq;
ok($seq->is_circular, 'Phi-X174 genome is circular');

# retrieving the spliced sequence from any split location requires
# spliced_seq()

my %sf_data = (
    #       start
    'A'  => [3981, 136, 1, 1542, 'join(3981..5386,1..136)', 'ATGGTTCGTT'],
    'A*' => [4497, 136, 1, 1026, 'join(4497..5386,1..136)', 'ATGAAATCGC'],
    'B'  => [5075, 136, 1, 363,  'join(5075..5386,1..51)',  'ATGGAACAAC'],
);

my @split_sfs = grep {
    $_->location->isa('Bio::Location::SplitLocationI')
    } $seq->get_SeqFeatures();

is(@split_sfs, 3, 'only 3 split locations');

for my $sf (@split_sfs) {
    isa_ok($sf->location, 'Bio::Location::SplitLocationI');
    my ($tag) = $sf->get_tag_values('product');
    my ($start, $end, $strand, $length, $ftstring, $first_ten) = @{$sf_data{$tag}};
    
    # these pass
    is($sf->location->to_FTstring, $ftstring, 'feature string');
    is($sf->spliced_seq->subseq(1,10), $first_ten, 'first ten nucleotides');
        is($sf->strand, $strand, 'strand');

    TODO: {
        local $TODO = "Need to define how to deal with start, end length for circular sequences";
        is($sf->start, $start, 'start');
        is($sf->end, $end, 'end');
        is($sf->length, $length, 'expected length');
    }
}

