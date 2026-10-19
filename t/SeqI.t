#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 4;

BEGIN {
    use_ok('Bio::SeqI');
}

{
    package MySeq;
    use base qw(Bio::SeqI);

    sub new {
        my $class = shift;
        return bless {}, $class;
    }

    sub display_id { return "test_id"; }
    sub seq        { return "ATGC"; }
    sub desc       { return "Test sequence"; }
    sub alphabet   { return "dna"; }
}

my $seq = MySeq->new();

isa_ok($seq, 'Bio::SeqI');
is($seq->display_id, 'test_id', 'display_id returns correct value');
is($seq->seq, 'ATGC', 'seq returns correct sequence');

