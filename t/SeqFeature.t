# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..27\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Computation;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

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

test 2, ( $feat->start == 40 );

test 3, ( $feat->end == 80 );

test 4, ( $feat->primary_tag eq 'exon' );


test 5, ( $feat->source_tag eq 'internal' );


$str = $feat->gff_string();
$str = ""; # shut up -w

# we need to figure out the correct mapping of this stuff
# soon

#if( $str ne "SEQ\tinternal\texon\t40\t80\t1\t.\t." ) {
#    print "not ok 3\n";
#} else {
#    print "ok 3\n";
#}

test 6, 1;

$pair = new Bio::SeqFeature::FeaturePair();

test 7, 1;

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

$pair->feature1($feat);
$pair->feature2($feat2);

test 8, 1;

test 9, ( $pair->start == 40 );

test 10, ( $pair->end == 80 );


test 11, ( $pair->primary_tag eq 'exon' );

test 12, ( $pair->source_tag eq 'internal' );

test 13, ( $pair->hstart == 400 );

test 14, ( $pair->hend == 440 );

test 15, ( $pair->hprimary_tag eq 'other' );

test 16, ( $pair->hsource_tag eq 'program_a' );

$pair->invert;

test 17, ( $pair->end == 440 );

# Test attaching a SeqFeature::Generic to a Bio::Seq
{
    # Make the parent sequence object
    my $seq = Bio::Seq->new(
        '-seq'          => 'aaaaggggtttt',
        '-display_id'   => 'test',
        '-moltype'      => 'dna',
        );
    
    # Make a SeqFeature
    my $sf1 = Bio::SeqFeature::Generic->new(
        '-start'    => 4,
        '-end'      => 9,
        '-strand'   => 1,
        );
    
    # Add the SeqFeature to the parent
    test 18, ($seq->add_SeqFeature($sf1));
    
    # Test that it gives the correct sequence
    my $sf_seq1 = $sf1->seq->seq;
    test 19, ($sf_seq1 eq 'aggggt');
    
    # Make a second seqfeature on the opposite strand
    my $sf2 = Bio::SeqFeature::Generic->new(
        '-start'    => 4,
        '-end'      => 9,
        '-strand'   => -1,
        );
    
    # This time add the PrimarySeq to the seqfeature
    # before adding it to the parent
    test 20, ($sf2->attach_seq($seq->primary_seq));
    $seq->add_SeqFeature($sf2);
    
    # Test again that we have the correct sequence
    my $sf_seq2 = $sf2->seq->seq;
    test 21, ($sf_seq2 eq 'acccct');
}

#Do some tests for computation.pm

test 22, ( $comp_obj1 = Bio::SeqFeature::Computation->new() );
test 23, ( $comp_obj1->computation_id(332) );
test 24, ( $comp_obj1->add_score_value('P', 33) );
{
    $comp_obj2 = Bio::SeqFeature::Computation->new();
    test 25, ($comp_obj1->add_sub_SeqFeature($comp_obj2, 'exon') );
    test 26, (@sft = $comp_obj1->all_sub_SeqFeature_types() );
    test 27, ($sft[0] eq 'exon');
}

