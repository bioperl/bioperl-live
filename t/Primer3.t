## $Id$

# test for Bio::Tools::Primer3.pm
# written by Rob Edwards

use strict;
use constant NUMTESTS => 10;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;

    plan tests => NUMTESTS;
}

use Bio::Tools::Primer3;
ok(1);

my ($p3, $num, $primer);

ok $p3=Bio::Tools::Primer3->new(-file=>"t/data/primer3_output.txt");
ok $num=$p3->number_of_results;
ok $num, 4, "Got $num";
ok $num=$p3->all_results;
ok defined $num, 1, "Can't get all results";
ok $num=$p3->primer_results(1);
ok defined $num, 1, "Can't get results for 1";
ok $primer=$p3->next_primer;
ok ref($primer) eq "Bio::Seq::PrimedSeq", 1, "reference for primer stream is not right";

