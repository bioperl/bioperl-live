# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 189);

    use_ok('Bio::Tools::IUPAC');
    use_ok('Bio::Seq');
    use_ok('Bio::PrimarySeq');
}


# IUPAC sequences and regular expressions

my $ambiseq = Bio::Seq->new(
    -seq      => 'ARTCGTTGN',
    -alphabet => 'dna',
);

my $ambiprimaryseq = Bio::Seq->new(
    -seq      => 'ARTCGTTGN',
    -alphabet => 'dna',
);

ok my $iupac = Bio::Tools::IUPAC->new( -seq => $ambiprimaryseq );

ok $iupac = Bio::Tools::IUPAC->new( -seq => $ambiseq );

ok my $regexp = $iupac->regexp;

is $regexp, 'A[AGR][TU]CG[TU][TU]G[ACGTUBDHKMNRSVWY]', 'Regexp';

is $iupac->count(), 80, 'Count';

my @seqs;
while (my $uniqueseq = $iupac->next_seq()) {
    push @seqs, $uniqueseq->seq;
    is $uniqueseq->isa('Bio::PrimarySeqI'), 1;
    like $uniqueseq->seq, qr/$regexp/i;
}

@seqs = sort @seqs;
is_deeply \@seqs, [
          'AATCGTTGA',
          'AATCGTTGC',
          'AATCGTTGG',
          'AATCGTTGT',
          'AATCGTTGU',
          'AATCGTUGA',
          'AATCGTUGC',
          'AATCGTUGG',
          'AATCGTUGT',
          'AATCGTUGU',
          'AATCGUTGA',
          'AATCGUTGC',
          'AATCGUTGG',
          'AATCGUTGT',
          'AATCGUTGU',
          'AATCGUUGA',
          'AATCGUUGC',
          'AATCGUUGG',
          'AATCGUUGT',
          'AATCGUUGU',
          'AAUCGTTGA',
          'AAUCGTTGC',
          'AAUCGTTGG',
          'AAUCGTTGT',
          'AAUCGTTGU',
          'AAUCGTUGA',
          'AAUCGTUGC',
          'AAUCGTUGG',
          'AAUCGTUGT',
          'AAUCGTUGU',
          'AAUCGUTGA',
          'AAUCGUTGC',
          'AAUCGUTGG',
          'AAUCGUTGT',
          'AAUCGUTGU',
          'AAUCGUUGA',
          'AAUCGUUGC',
          'AAUCGUUGG',
          'AAUCGUUGT',
          'AAUCGUUGU',
          'AGTCGTTGA',
          'AGTCGTTGC',
          'AGTCGTTGG',
          'AGTCGTTGT',
          'AGTCGTTGU',
          'AGTCGTUGA',
          'AGTCGTUGC',
          'AGTCGTUGG',
          'AGTCGTUGT',
          'AGTCGTUGU',
          'AGTCGUTGA',
          'AGTCGUTGC',
          'AGTCGUTGG',
          'AGTCGUTGT',
          'AGTCGUTGU',
          'AGTCGUUGA',
          'AGTCGUUGC',
          'AGTCGUUGG',
          'AGTCGUUGT',
          'AGTCGUUGU',
          'AGUCGTTGA',
          'AGUCGTTGC',
          'AGUCGTTGG',
          'AGUCGTTGT',
          'AGUCGTTGU',
          'AGUCGTUGA',
          'AGUCGTUGC',
          'AGUCGTUGG',
          'AGUCGTUGT',
          'AGUCGTUGU',
          'AGUCGUTGA',
          'AGUCGUTGC',
          'AGUCGUTGG',
          'AGUCGUTGT',
          'AGUCGUTGU',
          'AGUCGUUGA',
          'AGUCGUUGC',
          'AGUCGUUGG',
          'AGUCGUUGT',
          'AGUCGUUGU'
];

like $ambiseq->seq, qr/$regexp/i, 'Regexp matches ambiguous sequences';
like 'ARTCGTTGW', qr/$regexp/i;


# IUPAC code methods

my %iupac;
ok %iupac = $iupac->iupac_iub(), 'Nucleic IUPAC';
ok exists $iupac{'A'};
ok not exists $iupac{'Z'};

ok %iupac = $iupac->iupac_iub_amb();
ok exists $iupac{'N'};
ok not exists $iupac{'A'};

ok %iupac = $iupac->iupac_rev_iub();

ok %iupac = $iupac->iupac_iup(), 'Proteic IUPAC';
ok exists $iupac{'A'};
ok exists $iupac{'Z'};

ok %iupac = $iupac->iupac_iup_amb();
ok exists $iupac{'B'};
ok not exists $iupac{'A'};

ok %iupac = $iupac->iupac();
ok not(exists $iupac{'Z'});

ok %iupac = $iupac->iupac_amb();
ok not(exists $iupac{'A'});

ok %iupac = Bio::Tools::IUPAC->new->iupac_iup;

