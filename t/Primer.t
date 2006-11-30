## $Id$

# test for Bio::SeqFeature::Primer
# written by Rob Edwards

use strict;
use constant NUMTESTS => 18;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;

    plan tests => NUMTESTS;
}

use Bio::SeqFeature::Primer;
ok(1);

my ($primer, $location, $start, $end, $strand, $id, $tm, $tme);

ok $primer=Bio::SeqFeature::Primer->new(-seq=>'CTTTTCATTCTGACTGCAACG');
ok $primer->seq->seq eq "CTTTTCATTCTGACTGCAACG";
ok $primer->primary_tag eq "Primer";
ok $location=$primer->location(500);
ok $location==500;
ok $start=$primer->start(2);
ok $start == 2;
ok $end=$primer->end(19);
ok $end == 19;
ok $strand=$primer->strand(-1);
ok $strand == -1;
ok $id=$primer->display_id('test');
ok $id eq "test";
ok $tm = $primer->Tm;
ok $tme = $primer->Tm_estimate;
ok int($tm) == 52;
ok int($tme) == 58;
