# This is -*-Perl-*- code
# $Id$
use strict;

my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 8;
}


use Bio::Seq::LargeLocatableSeq;
use Data::Dumper;
ok 1;

ok my $llseq  = Bio::Seq::LargeLocatableSeq->new(-seq => 'at-cg',
                                                 -display_id => 'seq1');

print Dumper $llseq if $DEBUG;

ok $llseq->isa("Bio::Seq::LargeSeqI");

ok $llseq->seq, 'at-cg';
ok $llseq->add_sequence_as_string('atcc'), 9;

ok $llseq->start, 1;

ok $llseq->end, 8;
ok $llseq->length, 9;

1;
