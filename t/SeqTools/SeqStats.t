# -*-Perl-*- Test Harness script for Bioperl
# $Id: SeqStats.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 47);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::Tools::SeqStats');
}

my ($seqobj, $count, $seqobj_stats, $wt);

my $str = Bio::SeqIO->new(-file=> test_input_file('multifa.seq'), '-format' => 'Fasta');
$seqobj = $str->next_seq();
ok defined $seqobj, 'new Bio::Root::IO object';

my $seq_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
ok defined $seq_stats && $seq_stats, 'new Bio::Tools:SeqStats object';

# eg for DNA sequence
my $hash_ref = $seq_stats->count_monomers();  
is ( $hash_ref->{'A'}, 80 , 'count_monomers()');

$hash_ref = $seq_stats-> count_codons();  
ok defined $hash_ref && $hash_ref , 'count_codons()';

my $weight = $seq_stats->get_mol_wt();
ok defined $weight && $weight , 'get_mol_wt()' ;

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

# protein

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

# Issue 3185: https://redmine.open-bio.org/issues/3185

$seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT*',
                               -alphabet=>'protein', -id=>'test');
is($seqobj->seq, 'MQSERGITIDISLWKFETSKYYVT*');
$seqobj_stats  =  Bio::Tools::SeqStats->new('-seq' => $seqobj);
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is int $$wt[0], 2896;
is int $$wt[1], 2896;

# selenocysteine
ok $seqobj = Bio::PrimarySeq->new(-seq=>'MQSERGITIDISLWKFETSKYYVT',
                                  -alphabet=>'protein');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 2896 ;

# RNA

$seqobj = Bio::PrimarySeq->new(-seq=>'UYXUYNNYU', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 2768;
is &round($$wt[1]), 2891;

ok $seqobj = Bio::PrimarySeq->new(-seq=>'TGCCGTGTGTGCTGCTGCT', -alphabet=>'rna');
$wt = Bio::Tools::SeqStats->get_mol_wt($seqobj);
is &round($$wt[0]), 6104 ;

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

