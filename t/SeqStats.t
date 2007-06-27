# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 39);
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::SeqStats');
}

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG',
			       -alphabet=>'dna', -id=>'test');
ok $seqobj_stats = Bio::Tools::SeqStats->new(-seq=>$seqobj);
isa_ok $seqobj_stats, 'Bio::Tools::SeqStats';

$count = $seqobj_stats->count_monomers();  # for DNA sequence
is $count->{'A'}, 3;
is $count->{'C'}, 4;
is $count->{'G'}, 5;
is $count->{'T'}, 4;

$count = $seqobj_stats->count_codons();
is $count->{'ACT'}, 2;
is $count->{'GTG'}, 1;
is $count->{'GCG'}, 1;
is $count->{'TCA'}, 1;


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTACTTCA', -alphabet=>'dna',
			       -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
$wt = $seqobj_stats->get_mol_wt();  # for DNA sequence
is &round($$wt[0]), 2738 ;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACXACNNCA',
			       -alphabet=>'dna', -id=>'test');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 2693;
is &round($$wt[1]), 2813;


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG',
			       -alphabet=>'dna', -id=>'test');
$count = Bio::Tools::SeqStats->count_monomers($seqobj);  # for DNA sequence
is $count->{'A'}, 3;
is $count->{'C'}, 4;
is $count->{'G'}, 5;
is $count->{'T'}, 4;

$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT',
                               -alphabet=>'protein', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
$count = $seqobj_stats->count_monomers();  # for amino sequence
is $$count{'M'}, 1;
is $$count{'I'}, 3;
is $$count{'Y'}, 2;
is $$count{'T'}, 3;
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is int $$wt[0], 2896;
is int $$wt[1], 2896;

$seqobj = Bio::PrimarySeq->new(-seq=>'UYXUYNNYU', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 2768;
is &round($$wt[1]), 2891;

ok $seqobj = Bio::PrimarySeq->new(-seq=>'TGCCGTGTGTGCTGCTGCT', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 6104 ;

# selenocysteine
ok $seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT',
                                  -alphabet=>'protein');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 2896 ;

#
# hydropathicity aka "gravy" score
#

# normal seq (should succeed)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'MSFVLVAPDMLATAAADVVQIGSAVSAGS',
                                  -alphabet=>'protein');
my $gravy = Bio::Tools::SeqStats->hydropathicity($seqobj);
is int($gravy*1000), 1224;  # check to nearest 0.1%

# ambiguous sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'XXXB**BS', -alphabet=>'protein');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
like $@, qr/ambiguous amino acids/i;

# empty sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'', -alphabet=>'protein');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
like $@, qr/hydropathicity not defined/i;

# DNA sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'GATTACA', -alphabet=>'dna');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
like $@, qr/only meaningful for protein/;


#
# Extra functions
#

# perl does not have an explicit rounding function
sub round { return int ((shift @_) + 0.5 ) }
