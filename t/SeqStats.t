# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use Test;
use strict;

BEGIN { plan tests => 8 }

use Bio::PrimarySeq;
use Bio::SeqStats;

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG', -moltype=>'dna', -id=>'test');
$seqobj_stats  =  Bio::SeqStats->new($seqobj);

ok defined($seqobj_stats) && ref($seqobj_stats) && $seqobj_stats->isa('Bio::SeqStats');

$count = $seqobj_stats->count_monomers();  # for DNA sequence
ok $$count{'A'} == 3 && $$count{'C'} == 4 && $$count{'G'} == 5 && $$count{'T'} == 4  ;

$count = $seqobj_stats->count_codons();
ok $$count{'ACT'} == 2 && $$count{'GTG'} == 1 && $$count{'GCG'} == 1 && $$count{'TCA'} == 1 ;


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTACTTCA', -moltype=>'dna', -id=>'test');
$seqobj_stats  =  Bio::SeqStats->new($seqobj);
$wt = $seqobj_stats->get_mol_wt();  # for DNA sequence
ok $$wt[0] == 2976 ;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACXACNNCA', -moltype=>'dna', -id=>'test');
$wt = Bio::SeqStats->get_mol_wt($seqobj);
ok $$wt[0] == 2976 && $$wt[1] == 3099;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG', -moltype=>'dna', -id=>'test');
$count = Bio::SeqStats->count_monomers($seqobj);  # for DNA sequence
ok $$count{'A'} == 3 && $$count{'C'} == 4 && $$count{'G'} == 5 && $$count{'T'} == 4  ;

$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT', -moltype=>'protein', -id=>'test');
$seqobj_stats  =  Bio::SeqStats->new($seqobj);
$count = $seqobj_stats->count_monomers();  # for amino sequence
ok $$count{'M'} == 1 && $$count{'I'} == 3 && $$count{'Y'} == 2 && $$count{'T'} == 3 ;

$seqobj = Bio::PrimarySeq->new(-seq=>'UYXUYNNYU', -moltype=>'rna');
$wt = Bio::SeqStats->get_mol_wt($seqobj);
ok $$wt[0] == 3054 && $$wt[1] == 3177;

