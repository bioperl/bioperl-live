# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 5;
}
use Bio::Seq;
use Bio::SeqFeature::SimilarityPair;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Root::IO;

# test SimilarityPair

my $seq = (Bio::SeqIO->new('-format' => 'fasta',
			  '-file' => Bio::Root::IO->catfile("t","data","AAC12660.fa")))->next_seq();
ok defined( $seq) && $seq->isa('Bio::SeqI');
my $blast = Bio::SearchIO->new('-file'=>Bio::Root::IO->catfile("t","data","blast.report"), '-format' => 'blast');
ok defined ($blast) && $blast->isa('Bio::SearchIO');
my $r = $blast->next_result;
my $hit = $r->next_hit;
ok defined ($hit) && $hit->isa('Bio::Search::Hit::HitI'), 1, ' hit is ' . ref($hit);
my $sim_pair = $hit->next_hsp;
ok defined($sim_pair) && $sim_pair->isa('Bio::SeqFeatureI');
$seq->add_SeqFeature($sim_pair);
ok $seq->all_SeqFeatures() == 1;
