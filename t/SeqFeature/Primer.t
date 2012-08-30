# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 25);

    use_ok('Bio::SeqFeature::Primer');
    use_ok('Bio::PrimarySeq');
}

my ($primer, $location, $start, $end, $strand, $id, $tm, $tme);

# Initialize with a Bio::PrimarySeq
my $seq = Bio::PrimarySeq->new(-seq => 'CTTTTCATTCTGACTGCAACG');
ok $primer = Bio::SeqFeature::Primer->new(-sequence => $seq);
is $primer->seq->seq, 'CTTTTCATTCTGACTGCAACG';

# Initialize with a sequence string
ok $primer = Bio::SeqFeature::Primer->new(
    -seq => 'CTTTTCATTCTGACTGCAACG',
    -TARGET => '5,3',
);
is $primer->display_id, 'SeqFeature Primer object';
is $primer->seq->seq, 'CTTTTCATTCTGACTGCAACG';
is $primer->primary_tag, 'Primer';
ok $id = $primer->display_id('test');
is $id, 'test';
is $primer->{'-TARGET'}, '5,3';

# Coordinates
ok $location = $primer->location(500);
is $location, 500;
ok $start = $primer->start(2);
is $start, 2;
ok $end = $primer->end(19);
is $end, 19;
ok $strand = $primer->strand(-1);
is $strand, -1;

# Melting temperatures
ok $tm = $primer->Tm;
is int($tm), 52;
ok $tm = $primer->Tm(-salt => 0.05, -oligo => 0.0000001);
ok $tme = $primer->Tm_estimate;
is int($tme), 58;
ok $tm = $primer->Tm_estimate(-salt => 0.05);
