# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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
    plan tests => 17;
}
use Bio::PrimarySeq;

ok(1);

my $seq = Bio::PrimarySeq->new('-seq'              =>'TTGGTGGCGTCAACT',
			       '-display_id'       => 'new-id',
			       '-alphabet'          => 'dna',
			       '-accession_number' => 'X677667',
			       '-desc'             =>'Sample Bio::Seq object');
ok defined $seq;
ok $seq->isa('Bio::PrimarySeqI');
ok $seq->accession_number(), 'X677667';
ok $seq->seq(), 'TTGGTGGCGTCAACT';
ok $seq->display_id(), 'new-id';
ok $seq->alphabet(), 'dna';
ok $seq->is_circular(), undef;
ok $seq->is_circular(1);

my $trunc = $seq->trunc(1,4);
ok defined $trunc;
ok $trunc->length, 4;
ok $trunc->seq(), 'TTGG', "Expecting TTGG. Got ".$trunc->seq();

my $rev = $seq->revcom();
ok defined $rev; 

ok $rev->seq(), 'AGTTGACGCCACCAA', 'rev faile was ' . $rev->seq();

#
# Translate
#

my $aa = $seq->translate();
ok $aa->seq, 'LVAST', "Translation: ". $aa->seq;

$seq->seq('TTGGTGGCGTCAACTTAA');
$aa = $seq->translate(undef, undef, undef, undef, 1);

# tests for non-Methionin initiator codon (AGT) coding for M
ok $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# test for character '?' in the sequence string
ok $seq->seq('TTGGTGGCG?CAACT'), 'TTGGTGGCG?CAACT';
