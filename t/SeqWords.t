# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

## SeqWords.t, based on SeqStats.t
# Derek Gatherer, 11th November 2003

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { use lib 't'; }
    use Test;
    plan tests => 16;
}

use Bio::PrimarySeq;
use Bio::Tools::SeqWords;
use vars ('$DEBUG');

ok 1;

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTGACTGGC',
			       -alphabet=>'dna', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqWords->new(-seq=>$seqobj);

ok defined($seqobj_stats) && ref($seqobj_stats) &&
    $seqobj_stats->isa('Bio::Tools::SeqWords');

$count = $seqobj_stats->count_words(4);
ok $count->{'ACTG'}, 3;
ok $count->{'TGGC'}, 1;
ok $count->{'GTCA'}, 1;

$count = $seqobj_stats->count_overlap_words(4);
ok $count->{'ACTG'}, 3;
ok $count->{'TGGC'}, 2;
ok $count->{'GTCA'}, 1;
ok $count->{'GTGG'}, 1;

# now test a protein
$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVTIDISSLWKF',
                               -alphabet=>'protein', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqWords->new('-seq' => $seqobj);
ok defined($seqobj_stats) && ref($seqobj_stats) &&
    $seqobj_stats->isa('Bio::Tools::SeqWords');

$count = $seqobj_stats->count_words(4);
ok $count->{'MQSE'}, 1;
ok $count->{'LWKF'}, 1;
ok $count->{'IDIS'}, 2;

$count = $seqobj_stats->count_overlap_words(4);
ok $count->{'MQSE'}, 1;
ok $count->{'LWKF'}, 2;
ok $count->{'IDIS'}, 2;
