# -*-Perl-*- Test Harness script for Bioperl
# $Id: SeqWords.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 22);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Tools::SeqWords');
}

my ($seqobj, $count, $seqobj_stats, $wt);

my $str = Bio::SeqIO->new(-file=> test_input_file('multifa.seq'), '-format' => 'Fasta');
$seqobj= $str->next_seq();
ok defined $seqobj, 'new Bio::Root::IO object';

my $words = Bio::Tools::SeqWords->new('-seq' => $seqobj);
my $hash = $words->count_words(6);
ok (defined $words, 'new Bio::Tools::SeqWords object');
ok (defined $hash, 'count_words');

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
