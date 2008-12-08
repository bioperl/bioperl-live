# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 18);
	
    use_ok('Bio::SeqFeature::Primer');
}

my ($primer, $location, $start, $end, $strand, $id, $tm, $tme);

ok $primer=Bio::SeqFeature::Primer->new(-seq=>'CTTTTCATTCTGACTGCAACG');
is $primer->seq->seq, "CTTTTCATTCTGACTGCAACG";
is $primer->primary_tag, "Primer";
ok $location=$primer->location(500);
is $location,500;
ok $start=$primer->start(2);
is $start, 2;
ok $end=$primer->end(19);
is $end, 19;
ok $strand=$primer->strand(-1);
is $strand, -1;
ok $id=$primer->display_id('test');
is $id, "test";
ok $tm = $primer->Tm;
ok $tme = $primer->Tm_estimate;
is int($tm), 52;
is int($tme), 58;
