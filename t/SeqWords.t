# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 19);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::SeqWords');
}

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTGACTGGC',
			       -alphabet=>'dna', -id=>'test');
ok $seqobj_stats  =  Bio::Tools::SeqWords->new(-seq=>$seqobj);
isa_ok $seqobj_stats, 'Bio::Tools::SeqWords';

$count = $seqobj_stats->count_words(4);
is $count->{'ACTG'}, 3;
is $count->{'TGGC'}, 1;
is $count->{'GTCA'}, 1;

$count = $seqobj_stats->count_overlap_words(4);
is $count->{'ACTG'}, 3;
is $count->{'TGGC'}, 2;
is $count->{'GTCA'}, 1;
is $count->{'GTGG'}, 1;

# now test a protein
$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVTIDISSLWKF',
                               -alphabet=>'protein', -id=>'test');
ok $seqobj_stats  =  Bio::Tools::SeqWords->new('-seq' => $seqobj);
isa_ok $seqobj_stats, 'Bio::Tools::SeqWords';

$count = $seqobj_stats->count_words(4);
is $count->{'MQSE'}, 1;
is $count->{'LWKF'}, 1;
is $count->{'IDIS'}, 2;

$count = $seqobj_stats->count_overlap_words(4);
is $count->{'MQSE'}, 1;
is $count->{'LWKF'}, 2;
is $count->{'IDIS'}, 2;
