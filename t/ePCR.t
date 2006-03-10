# -*-Perl-*-
## Bioperl Test Harness Script for Modules

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
    use vars qw($NUMTESTS $DEBUG);
    $NUMTESTS = 25;
    $DEBUG   = 1;
    plan test => $NUMTESTS;
}

use Bio::Tools::EPCR;
use Bio::SeqIO;
use Bio::Root::IO;

my $seqio = new Bio::SeqIO('-format' => 'fasta', '-file' => Bio::Root::IO->catfile("t","data","genomic-seq.fasta"));

my $seq = $seqio->next_seq;
ok($seq);
my $epcr = new Bio::Tools::EPCR( '-file' => Bio::Root::IO->catfile("t","data","genomic-seq.epcr"));
ok ($epcr);
my %strand;
while( defined(my $feature = $epcr->next_feature) ) {
    ok($feature);
    ok($feature->start);
    ok($feature->end);
    $seq->add_SeqFeature($feature);
    $strand{$feature->strand} ++;
}
ok ($strand{1},  3 , 'expected 3 forward strand ePCR hits');
ok ($strand{-1}, 3 , 'expected 3 reverse strand ePCR hits');


if( $DEBUG ) {
    $seqio = new Bio::SeqIO('-format' => 'genbank' );
    $seqio->write_seq($seq);
}
