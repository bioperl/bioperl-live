# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 33);

    use_ok('Bio::SeqFeature::Primer');
    use_ok('Bio::PrimarySeq');
}

my ($primer, $location, $start, $end, $strand, $id, $tm, $tme);

# Initialize with a Bio::PrimarySeq
my $seq = Bio::PrimarySeq->new(-seq => 'CTTTTCATTCTGACTGCAACG');
ok $primer = Bio::SeqFeature::Primer->new(-sequence => $seq);
is $primer->seq->isa('Bio::PrimarySeqI'), 1;
is $primer->seq->seq, 'CTTTTCATTCTGACTGCAACG';

# Initialize with a sequence string
ok $primer = Bio::SeqFeature::Primer->new(
    -seq    => 'CTTTTCATTCTGACTGCAACG',
    -TARGET => '5,3',
    -start  => 3,
);
is $primer->start, 3;
is $primer->display_name, 'SeqFeature Primer object';
is $primer->seq->isa('Bio::PrimarySeqI'), 1;
is $primer->seq->seq, 'CTTTTCATTCTGACTGCAACG';
is $primer->primary_tag, 'Primer';
ok $id = $primer->display_name('test');
is $id, 'test';
is $primer->{'-TARGET'}, '5,3';

# Coordinates
ok $start = $primer->start(2);
is $start, 2;
ok $end = $primer->end(19);
is $end, 19;
ok $strand = $primer->strand(-1);
is $strand, -1;

# Legacy: passing a string to location
{
   local $SIG{'__WARN__'} = sub {  }; # Silence all warnings (we expect deprecation messages)
   ok $location = $primer->location(500);
   isa_ok $location, 'Bio::Location::Simple';
   is $primer->start, 500;
   ok $location = $primer->location('3,25');
   isa_ok $location, 'Bio::Location::Simple';
   is $primer->start, 3;
   is $primer->end, 25;
}

# Melting temperatures
ok $tm = $primer->Tm;
is int($tm), 52;
ok $tm = $primer->Tm(-salt => 0.05, -oligo => 0.0000001);
ok $tme = $primer->Tm_estimate;
is int($tme), 58;
ok $tm = $primer->Tm_estimate(-salt => 0.05);


#####
use Data::Dumper;
print Dumper($primer);
#####
