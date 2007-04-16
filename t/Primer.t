## $Id$

# test for Bio::SeqFeature::Primer
# written by Rob Edwards

use strict;
use constant NUMTESTS => 18;

BEGIN {
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;

    plan tests => NUMTESTS;
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
