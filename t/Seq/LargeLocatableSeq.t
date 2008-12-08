# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 8);
	
	use_ok('Bio::Seq::LargeLocatableSeq');
}

ok my $llseq  = Bio::Seq::LargeLocatableSeq->new(-seq => 'at-cg',
                                                 -display_id => 'seq1');

isa_ok $llseq, "Bio::Seq::LargeSeqI";

is $llseq->seq, 'at-cg';
is $llseq->add_sequence_as_string('atcc'), 9;

is $llseq->start, 1;

is $llseq->end, 8;
is $llseq->length, 9;
