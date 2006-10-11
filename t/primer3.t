# -*-Perl-*- mode (to keep my emacs happy)
## $Id$

# test for Bio::Tools::Primer3.pm
# written by Rob Edwards
# and Chad Matsalla

use strict;
use vars qw($NUMTESTS $DEBUG $ERROR $XML_ERROR);


BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;
    $NUMTESTS  = 24;

    plan tests => $NUMTESTS;

    eval {  require Clone; };
    if ( $@ ) {
	warn("Clone not installed. This means that the module is not usable. Skipping tests\n");
	$ERROR = 1;
    }
}

END {
        foreach ( $Test::ntest..$NUMTESTS) {
	skip("Missing dependencies. Skipping tests",1);
    }
}

exit 0 if $ERROR;

require Bio::Tools::Primer3;
ok(1);

my ($p3, $num, $primer);

ok $p3 = Bio::Tools::Primer3->new(-file => File::Spec->catfile(qw(t data primer3_output.txt)));
ok $num = $p3->number_of_results;
ok $num, 5, "Got $num";
ok $num = $p3->all_results;
ok defined $num, 1, "Can't get all results";
ok $num = $p3->primer_results(1);
ok defined $num, 1, "Can't get results for 1";
ok $primer = $p3->next_primer;
ok ref($primer) eq "Bio::Seq::PrimedSeq", 1, 
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
ok $left_primer_seq->seq, "GAGGGTAACACGCTGGTCAT";
