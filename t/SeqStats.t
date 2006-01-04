# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use strict;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { use lib 't'; }
    use Test;
    plan tests => 36;
}

use Bio::PrimarySeq;
use Bio::Tools::SeqStats;
use vars ('$DEBUG');

my ($seqobj, $count, $seqobj_stats, $wt);

$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG',
			       -alphabet=>'dna', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new(-seq=>$seqobj);

ok defined($seqobj_stats) && ref($seqobj_stats) &&
    $seqobj_stats->isa('Bio::Tools::SeqStats');

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


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTACTTCA', -alphabet=>'dna',
			       -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
$wt = $seqobj_stats->get_mol_wt();  # for DNA sequence
ok &round($$wt[0]), 2738 ;

$seqobj = Bio::PrimarySeq->new(-seq=>'ACXACNNCA',
			       -alphabet=>'dna', -id=>'test');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok &round($$wt[0]), 2693;
ok &round($$wt[1]), 2813;


$seqobj = Bio::PrimarySeq->new(-seq=>'ACTGTGGCGTCAACTG',
			       -alphabet=>'dna', -id=>'test');
$count = Bio::Tools::SeqStats->count_monomers($seqobj);  # for DNA sequence
ok $count->{'A'}, 3;
ok $count->{'C'}, 4;
ok $count->{'G'}, 5;
ok $count->{'T'}, 4;

$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT',
                               -alphabet=>'protein', -id=>'test');
$seqobj_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
$count = $seqobj_stats->count_monomers();  # for amino sequence
ok $$count{'M'}, 1;
ok $$count{'I'}, 3;
ok $$count{'Y'}, 2;
ok $$count{'T'}, 3;
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok int $$wt[0], 2896;
ok int $$wt[1], 2896;

$seqobj = Bio::PrimarySeq->new(-seq=>'UYXUYNNYU', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok &round($$wt[0]), 2768;
ok &round($$wt[1]), 2891;

ok $seqobj = Bio::PrimarySeq->new(-seq=>'TGCCGTGTGTGCTGCTGCT', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok &round($$wt[0]), 6104 ;

# selenocysteine
ok $seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT',
                                  -alphabet=>'protein');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
ok &round($$wt[0]), 2896 ;

#
# hydropathicity aka "gravy" score
#

# normal seq (should succeed)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'MSFVLVAPDMLATAAADVVQIGSAVSAGS',
                                  -alphabet=>'protein');
my $gravy = Bio::Tools::SeqStats->hydropathicity($seqobj);
ok int($gravy*1000) == 1224;  # check to nearest 0.1%

# ambiguous sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'XXXB**BS', -alphabet=>'protein');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
ok $@ =~ /ambiguous amino acids/i;

# empty sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'', -alphabet=>'protein');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
ok $@ =~ /hydropathicity not defined/i;

# DNA sequence (should fail)
ok $seqobj = Bio::PrimarySeq->new(-seq=>'GATTACA', -alphabet=>'dna');
eval { Bio::Tools::SeqStats->hydropathicity($seqobj) };
ok $@ =~ /only meaningful for protein/;


#
# Extra functions
#

# perl does not have an explicit rounding function
sub round { return int ((shift @_) + 0.5 ) }
