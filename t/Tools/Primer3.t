# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 14);

    use_ok('Bio::Tools::Primer3');
}

my ($p3, $num, $primer);

ok $p3 = Bio::Tools::Primer3->new(-file => test_input_file('primer3_output.txt'));
ok $num = $p3->number_of_results;
is $num, 5 or diag "Got $num";
ok $num = $p3->all_results;
is defined $num, 1 or diag "Can't get all results";
ok $num = $p3->primer_results(1);
is defined $num, 1 or diag "Can't get results for 1";
ok $primer = $p3->next_primer;
isa_ok $primer, "Bio::Seq::PrimedSeq" or diag
  "reference for primer stream is not right";

# get the left primer
my $left_primer = $primer->get_primer('left');

# get the sequence for that primer. This is a test to verify behavior 
# on the bioperl list in or about 050315
my $seqobj = $left_primer->seq();

my $seq = $seqobj->seq();

my $other_left_primer = $primer->get_primer();

# a different way to access the primers in the stream
my $alt = $p3->primer_results(0,'PRIMER_LEFT_INPUT');

# next one
ok $primer = $p3->next_primer;
# get the left primer
my $left_primer_seq = $primer->get_primer('left')->seq;
is $left_primer_seq->seq, "GAGGGTAACACGCTGGTCAT";

# bug 2862
ok $p3 = Bio::Tools::Primer3->new(-file => test_input_file('bug2862.pmr'));
$num = 0;
while ($p3->next_primer) { $num++ };
is $p3->number_of_results, $num, 'bug 2862';
