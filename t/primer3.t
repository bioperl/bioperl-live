## $Id$

# test for Bio::Tools::Primer3.pm
# written by Rob Edwards
# and Chad Matsalla

use strict;
use Dumpvalue();
my $dumper = new Dumpvalue();

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
     # print("This is the primer object: ($primer)\n");
     # get the left primer
my $left_primer = $primer->get_primer('left');
     # print("This is the left primer object: ($left_primer)\n");
     # get the sequence for that primer. This is a test to verify behavior on the bioperl list in or about 050315
my $seqobj = $left_primer->seq();
     # print("This is the sequence object for the left primer: ($seqobj)\n");
my $seq = $seqobj->seq();
     # print("This is the sequence for the left primer: ($seq)\n");
my $other_left_primer = $primer->get_primer();
     # a different way to access the primers in the stream
     # print("This is the 0th primer result:\n");
my $alt = $p3->primer_results(0,'PRIMER_LEFT_INPUT');
     # $dumper->dumpValue($alt);

