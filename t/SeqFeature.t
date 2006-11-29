# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS);
my $skipdbtests;
my $skip_all;
BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	$NUMTESTS = 211;
	plan tests => $NUMTESTS;

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

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('Skipping tests which need the Bio::DB::GenBank module',1);
	}
}

exit(0) if $skip_all;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::SimilarityPair;
use Bio::SeqFeature::Computation;
require Bio::SeqFeature::Annotated;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::UTR;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Poly_A_site;
use Bio::SeqFeature::Gene::GeneStructure;

use Bio::Location::Fuzzy;
use Env qw(BIOPERLDEUG); # for importing bioperldebug var
ok(1);

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

ok $feat->start, 40;

ok $feat->end, 80;

ok $feat->primary_tag, 'exon';

ok $feat->source_tag, 'internal';

$str = $feat->gff_string() || ""; # placate -w

# we need to figure out the correct mapping of this stuff
# soon

#if( $str ne "SEQ\tinternal\texon\t40\t80\t1\t.\t." ) {
#    print "not ok 3\n";
#} else {
#    print "ok 3\n";
#}

ok(1);

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

ok $pair->feature1, $feat;
ok $pair->feature2, $feat2;
ok $pair->start, 40;
ok $pair->end, 80;
ok $pair->primary_tag, 'exon';
ok $pair->source_tag, 'internal';
ok $pair->hstart, 400;
ok $pair->hend, 440;
ok $pair->hprimary_tag, 'other' ;
ok $pair->hsource_tag, 'program_a';

$pair->invert;
ok $pair->end, 440;

# Test attaching a SeqFeature::Generic to a Bio::Seq
{
    # Make the parent sequence object
    my $seq = Bio::Seq->new(
        '-seq'          => 'aaaaggggtttt',
        '-display_id'   => 'test',
        '-alphabet'      => 'dna',
        );
    
    # Make a SeqFeature
    my $sf1 = Bio::SeqFeature::Generic->new(
        '-start'    => 4,
        '-end'      => 9,
        '-strand'   => 1,
        );
    
    # Add the SeqFeature to the parent
    ok ($seq->add_SeqFeature($sf1));
    
    # Test that it gives the correct sequence
    my $sf_seq1 = $sf1->seq->seq;
    ok $sf_seq1, 'aggggt';
    ok $sf1->end,9;
    ok $sf1->start,4;

    # Make a second seqfeature on the opposite strand
    my $sf2 = Bio::SeqFeature::Generic->new(
        '-start'    => 4,
        '-end'      => 9,
        '-strand'   => -1,
        );
    
    # This time add the PrimarySeq to the seqfeature
    # before adding it to the parent
    ok ($sf2->attach_seq($seq->primary_seq));
    $seq->add_SeqFeature($sf2);
    
    # Test again that we have the correct sequence
    my $sf_seq2 = $sf2->seq->seq;
    ok $sf_seq2, 'acccct';
}

#Do some tests for computation.pm

ok defined ( $comp_obj1 = Bio::SeqFeature::Computation->new('-start' => 1,
							    '-end'   => 10) );
ok ( $comp_obj1->computation_id(332),332 );
ok ( $comp_obj1->add_score_value('P', 33) );
{
    $comp_obj2 = Bio::SeqFeature::Computation->new('-start' => 2,
						   '-end'   => 10);
    ok ($comp_obj1->add_sub_SeqFeature($comp_obj2, 'exon') );
    ok (@sft = $comp_obj1->all_sub_SeqFeature_types() );
    ok ($sft[0], 'exon');
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

ok ( $comp_obj1->computation_id, 2231 );
ok ( $comp_obj1->add_score_value('P', 33) );
ok ( ($comp_obj1->each_score_value('no_score'))[0], '334');

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

ok($sfeat->start, 2);
ok($sfeat->end, 8);

# some tests to see if we can set a feature to start at 0
$sfeat = new Bio::SeqFeature::Generic(-start => 0, -end => 0 );

ok(defined $sfeat->start);
ok($sfeat->start,0);
ok(defined $sfeat->end);
ok($sfeat->end,0);


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
ok($t->mrna->length, 1693);
ok($gene->utrs, 2);



# Test for bug when Locations are not created explicitly

my $feat1 = new Bio::SeqFeature::Generic(-start => 1,
					 -end   => 15,
					 -strand=> 1);

$feat2 = new Bio::SeqFeature::Generic(-start => 10,
					 -end   => 25,
					 -strand=> 1);

my $overlap = $feat1->location->union($feat2->location);
ok($overlap->start, 1);
ok($overlap->end,   25);

my $intersect = $feat1->location->intersection($feat2->location);
ok($intersect->start, 10);
ok($intersect->end,   15);


# now let's test spliced_seq

ok  $seqio = new Bio::SeqIO(-file => Bio::Root::IO->catfile
			    (qw(t data AY095303S1.gbk)),
                            -format  => 'genbank');

ok $geneseq = $seqio->next_seq();
my ($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
my $db;

unless( $skipdbtests ) {
 $db = new Bio::DB::GenBank(-verbose=> $ENV{BIOPERLDEBUG});
 $CDS->verbose(-1);
 my $cdsseq = $CDS->spliced_seq(-db => $db,-nosort => 1);
 
 ok($cdsseq->subseq(1,60, 'ATGCAGCCATACGCTTCCGTGAGCGGGCGATGTCTATC'.
		    'TAGACCAGATGCATTGCATGTGATACCGTTTGGGCGAC'));
 ok($cdsseq->translate->subseq(1,100), 'MQPYASVSGRCLSRPDALHVIPFGRP'.
    'LQAIAGRRFVRCFAKGGQPGDKKKLNVTDKLRLGNTPPTLDVLKAPRPTDAPSAIDDAPSTSGLGLGGGVASPR');
} else {
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',1);
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',1);

}
ok  $seqio = new Bio::SeqIO(-file => Bio::Root::IO->catfile
			    (qw(t data AF032047.gbk)),
                            -format  => 'genbank');
ok $geneseq = $seqio->next_seq();
($CDS) = grep { $_->primary_tag eq 'CDS' } $geneseq->get_SeqFeatures;
unless ($skipdbtests ) {
    my $cdsseq = $CDS->spliced_seq( -db => $db, -nosort => 1);
    ok($cdsseq->subseq(1,60, 'ATGGCTCGCTTCGTGGTGGTAGCCCTGCTCGCGCTACTCTCTCTG'.
		       'TCTGGCCTGGAGGCTATCCAGCATG'));
    ok($cdsseq->translate->seq, 'MARFVVVALLALLSLSGLEAIQHAPKIQVYSRHPAENGKPNFL'.
       'NCYVSGFHPSDIEVDLLKNGKKIEKVEHSDLSFSKDWSFYLLYYTEFTPNEKDEYACRVSHVTFPTPKTVKWDRTM*');
} else {
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',1);
    skip('Cannot test for remote loc with spliced_seq w/o LWP installed',1);
}


# trans-spliced 

ok( $seqio = Bio::SeqIO->new(-format => 'genbank',
									  -file   => 
			    Bio::Root::IO->catfile(qw(t data NC_001284.gbk))));
my $genome = $seqio->next_seq;

foreach my $cds (grep { $_->primary_tag eq 'CDS' } $genome->get_SeqFeatures) {
   my $spliced = $cds->spliced_seq(-nosort => 1)->translate->seq;
   chop($spliced); # remove stop codon
   ok($spliced,($cds->get_tag_values('translation'))[0],'spliced seq translation matches expected');
}

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
{
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
	eval {
		$sfa2 = Bio::SeqFeature::Annotated->new(-feature => $sfg);
	};
	if ($@) {
		foreach ( $Test::ntest..$NUMTESTS ) { skip('Could not get sofa definitions from external server',1); }
		exit(0);
	}
  ok $sfa2->type->name,'nucleotide_motif';
  ok $sfa2->primary_tag, 'nucleotide_motif';
  ok $sfa2->source,'program_a';
  ok $sfa2->strand,1;
  ok $sfa2->start,400;
  ok $sfa2->end,440;
  ok $sfa2->get_Annotations('silly')->value,20;
  ok $sfa2->get_Annotations('new')->value,1;

  my $sfa3 = Bio::SeqFeature::Annotated->new( -start => 1,
					      -end => 5,
					      -strand => "+",
					      -frame => 2,
					      -phase => 2,
					      -score => 12,
					      -display_name => 'test.annot',
					      -seq_id => 'test.displayname' );
  eval {
	$sfa3->from_feature($sfg);
  };
  if ($@) {
	foreach ( $Test::ntest..$NUMTESTS ) { skip('Could not get sofa definitions from external server',1); }
	exit(0);
  }
  ok $sfa3->type->name,'nucleotide_motif';
  ok $sfa3->primary_tag, 'nucleotide_motif';
  ok $sfa3->source,'program_a';
  ok $sfa3->strand,1;
  ok $sfa3->start,400;
  ok $sfa3->end,440;
  ok $sfa3->get_Annotations('silly')->value,20;
  ok $sfa3->get_Annotations('new')->value,1;
}
