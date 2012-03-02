# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 39);
    use_ok('Bio::SeqFeature::Primer');
    use_ok('Bio::PrimarySeq');
}

my ($primer, $primer_seq, $location, $start, $end, $strand, $id, $tm, $tme, $template, $seq);

# Implied primer sequence
$template = Bio::Seq->new( -seq => 'AAAAACCCCCGGGGGTTTTT' );
ok $primer = Bio::SeqFeature::Primer->new( -start => 6, -end => 10 ), 'Implied primer sequence';
ok $primer->attach_seq($template); # attach the primer to a template
ok $primer_seq = $primer->seq;
isa_ok $primer_seq, 'Bio::PrimarySeqI';
is $primer_seq->seq, 'CCCCC';

# Bio::PrimarySeq primer
$template = Bio::Seq->new( -seq => 'AAAAACCCCCGGGGGTTTTT' );
$seq = Bio::PrimarySeq->new(-seq => 'CTTTTCATTCTGACTGCAACG');
ok $primer = Bio::SeqFeature::Primer->new(-seq => $seq), 'PrimarySeq primer';
ok $primer->attach_seq($template);
ok $primer_seq = $primer->seq;
isa_ok $primer_seq, 'Bio::PrimarySeqI';
is $primer_seq->seq, 'CTTTTCATTCTGACTGCAACG';

# Initialize with a sequence string
ok $primer = Bio::SeqFeature::Primer->new(
    -seq    => 'CTTTTCATTCTGACTGCAACG',
    -start  => 3,
    -id     => 'myid',
);
is $primer->start, 3;
ok $primer_seq = $primer->seq;
is $primer_seq->isa('Bio::PrimarySeqI'), 1;
is $primer_seq->seq, 'CTTTTCATTCTGACTGCAACG';
is $primer_seq->id, 'myid';
is $primer->primary_tag, 'Primer';
ok $primer->display_name('test');
is $primer->display_name, 'test';

# Coordinates
ok $primer->start(2);
is $primer->start, 2;
ok $primer->end(19);
is $primer->end, 19;
ok $primer->strand(-1);
is $primer->strand, -1;
ok $location = $primer->location;
isa_ok $location, 'Bio::LocationI';

# Melting temperatures
ok $tm = $primer->Tm;
is int($tm), 52;
ok $tm = $primer->Tm(-salt => 0.05, -oligo => 0.0000001);
ok $tme = $primer->Tm_estimate;
is int($tme), 58;
ok $tm = $primer->Tm_estimate(-salt => 0.05);

# Legacy
#   * initializing with -sequence
#   * passing a string to location
{
   local $SIG{'__WARN__'} = sub {  }; # Silence deprecation warnings
   ok $primer = Bio::SeqFeature::Primer->new(
       -sequence => 'CTTTTCATTCTGACTGCAACG',
       -TARGET   => '5,3',
       -start    => 3,
   );
   is $primer->{'-TARGET'}, '5,3';
   ok $location = $primer->location('3,25');
   is $location, '3,25';
}





