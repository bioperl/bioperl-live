# -*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    use vars qw($NUMTESTS $DEBUG);
    $NUMTESTS = 28;
    $DEBUG   = $ENV{'BIOPERLDEBUG'} || 0;
    plan tests => $NUMTESTS;
    use_ok('Bio::Tools::EPCR');
    use_ok('Bio::SeqIO');
    use_ok('Bio::Root::IO');
}

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
is ($strand{1},  3, 'got 3 forward strand ePCR hits');
is ($strand{-1}, 3, 'got 3 reverse strand ePCR hits');

if( $DEBUG ) {
    $seqio = new Bio::SeqIO('-format' => 'genbank' );
    $seqio->write_seq($seq);
}
