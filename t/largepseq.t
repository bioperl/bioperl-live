# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 7 }
use Bio::Seq::LargePrimarySeq;

my $pseq = Bio::Seq::LargePrimarySeq->new();
ok $pseq;
$pseq->add_sequence_as_string('ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAAT');
$pseq->add_sequence_as_string('GTTTGGGGTTAAACCCCTTTGGGGGGT');

ok $pseq->display_id('hello'), 'hello';

ok $pseq->seq, 'ATGGGGTGGGGTGAAACCCTTTGGGGGTGGGGTAAATGTTTGGGGTTAAACCCCTTTGGGGGGT' , "Sequence is " . $pseq->seq;

ok $pseq->subseq(3,7), 'GGGGT', "Subseq is ".$pseq->subseq(3,7);
ok ($pseq->trunc(8,15)->seq, 'GGGGTGAA', 
    'trunc seq was ' . $pseq->trunc(8,15)->seq);


ok $pseq->moltype('dna'), 'dna'; # so translate will not complain
ok $pseq->translate()->seq, 'MGWGETLWGWGKCLGLNPFGG';
