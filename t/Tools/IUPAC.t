# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 42);

    use_ok('Bio::Tools::IUPAC');
    use_ok('Bio::Seq');
}


# IUPAC sequences and regular expressions

my $ambiseq = Bio::Seq->new(
    -seq      => 'ARTCGTTGN',
    -alphabet => 'dna',
);

my $iupac = Bio::Tools::IUPAC->new( -seq => $ambiseq );

ok my $regexp = $iupac->regexp;

like $regexp, qr/A\[AGR\]TCGTTG\[ACGTBDHKMNRSVWY\]/, 'Regexp';

is $iupac->count(), 8, 'Count';

my @seqs;
while (my $uniqueseq = $iupac->next_seq()) {
    push @seqs, $uniqueseq->seq;
    is $uniqueseq->isa('Bio::PrimarySeqI'), 1;
    like $uniqueseq->seq, $regexp;
}

@seqs = sort @seqs;
is_deeply \@seqs, [ 'AATCGTTGA', 'AATCGTTGC', 'AATCGTTGG', 'AATCGTTGT',
                    'AGTCGTTGA', 'AGTCGTTGC', 'AGTCGTTGG', 'AGTCGTTGT' ];

like $ambiseq->seq, $regexp, 'Regexp matches ambiguous sequences';
like 'ARTCGTTGW', $regexp;


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

