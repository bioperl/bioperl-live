# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);

    use_ok('Bio::Tools::IUPAC');
    use_ok('Bio::Seq');
}

# test IUPAC

my $ambiseq = Bio::Seq->new(-seq => 'ARTCGTTGN',
                            -alphabet => 'dna'); 

my $stream = Bio::Tools::IUPAC->new(-seq => $ambiseq);
is $stream->count(), 8;

my @seqs;
my $valid_obj = 1;
while (my $uniqueseq = $stream->next_seq()) {
    push @seqs, $uniqueseq->seq;
    if( ! $uniqueseq->isa('Bio::Seq') ) {
        $valid_obj = 0;
        last; # no point continuing if we get here
    }
}
ok $valid_obj, 'Valid Bio::Seq objects';

@seqs = sort @seqs;
is_deeply \@seqs, [ 'AATCGTTGA', 'AATCGTTGC', 'AATCGTTGG', 'AATCGTTGT',
                    'AGTCGTTGA', 'AGTCGTTGC', 'AGTCGTTGG', 'AGTCGTTGT' ];
