# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use Test;
use strict;

BEGIN { plan tests => 22}

use Bio::PrimarySeq;
use Bio::Tools::SeqStats;

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG', -moltype=>'dna', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new($seqobj);

ok defined($seqobj_stats) && ref($seqobj_stats) && $seqobj_stats->isa('Bio::Tools::SeqStats');

$count = $seqobj_stats->count_monomers();  # for DNA sequence
ok $count->{'A'}, 3; 
ok $count->{'C'}, 4; 
ok $count->{'G'}, 5;
ok $count->{'T'}, 4;

$count = $seqobj_stats->count_codons();
ok $count->{'ACT'}, 2;
ok $count->{'GTG'}, 1;
ok $count->{'GCG'}, 1;
ok $count->{'TCA'}, 1;


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTACTTCA', -moltype=>'dna', 
			       -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new($seqobj);
$wt = $seqobj_stats->get_mol_wt();  # for DNA sequence
ok $wt->[0], 2976;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACXACNNCA', 
			       -moltype=>'dna', -id=>'test');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok $wt->[0], 2976; 
ok $wt->[1], 3099;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG', 
			       -moltype=>'dna', -id=>'test');
$count = Bio::Tools::SeqStats->count_monomers($seqobj);  # for DNA sequence
ok $count->{'A'}, 3;
ok $count->{'C'}, 4; 
ok $count->{'G'}, 5;
ok $count->{'T'}, 4;

$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT', -moltype=>'protein', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new($seqobj);
$count = $seqobj_stats->count_monomers();  # for amino sequence
ok $count->{'M'}, 1;
ok $count->{'I'}, 3; 
ok $count->{'Y'}, 2; 
ok $count->{'T'}, 3;

$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok $wt->[0], 2896;
ok $wt->[0], 2896;

$seqobj = Bio::PrimarySeq->new(-seq=>'UYXUYNNYU', -moltype=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok $wt->[0], 3054;
ok $wt->[0], 3177;

