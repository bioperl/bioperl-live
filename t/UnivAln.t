# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 11;
}
use Bio::UnivAln;

ok(1);
my($s,@s);

my $aln = Bio::UnivAln->new(-seqs=>"TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTG--TCAACTG");
ok $aln;
ok $aln->layout("fasta");
ok( $aln->seqs(1,1),  "TCCCGCGTCAACTG\n");
ok( $aln->seqs(2),  "TGGTGCTTCAACCG\nACTTG--TCAACTG\n");
ok( $aln->seqs([2..3]), "TGGTGCTTCAACCG\nACTTG--TCAACTG\n");
ok( $aln->consensus(), '!!!!G!!TCAAC!G' );
ok( $aln->complement(1), '' );
ok( $aln->gap_free_sites(), "TCCCGTCAACTG\nTGGTGTCAACCG\nACTTGTCAACTG\n");
ok( $aln->no_allgap_sites(), "TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTG--TCAACTG\n");
ok( $aln->remove_gaps(), "TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTGTCAACTG\n")
