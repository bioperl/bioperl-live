# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 46);

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

ok my $regexp = $iupac->regexp, 'Regexp';
is $regexp, 'A[AGR]TCGTTG[ACGTBDHKMNRSVWY]';

$regexp = $iupac->regexp(1);
is $regexp, 'A[AGR][TU]CG[TU][TU]G[ACGTUBDHKMNRSVWY]', 'Regexp';


is $iupac->count(), 8, 'Count';

my @seqs;
while (my $uniqueseq = $iupac->next_seq()) {
    push @seqs, $uniqueseq->seq;
    is $uniqueseq->isa('Bio::PrimarySeqI'), 1;
    like $uniqueseq->seq, qr/$regexp/i;
}

@seqs = sort @seqs;
is_deeply \@seqs, [ 'AATCGTTGA', 'AATCGTTGC', 'AATCGTTGG', 'AATCGTTGT',
                    'AGTCGTTGA', 'AGTCGTTGC', 'AGTCGTTGG', 'AGTCGTTGT' ];

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

