# -*-Perl-*-
## Bioperl Test Harness Script for various modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..4\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::Tools::Sigcleave;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

sub test ($$;$) {
    my($num, $true,$msg) = @_;
    print($true ? "ok $num\n" : "not ok $num $msg\n");
}

$protein = 'MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN';

## need this to escape the strict() warnings...
$formatted_output = "";
my %results;

# Build object
test 2, $sigcleave_object = new Bio::Tools::Sigcleave(-id         =>'test_sigcleave_seq',
                                                      -type      =>'amino',
                                                      -threshold => 0,
			                              -seq       =>$protein); 
# Test raw result accessor
test 3, %results = $sigcleave_object->signals, "unable to get raw sigcleave results";

# Test formatted output method
test 4, $formatted_output = $sigcleave_object->pretty_print, "unable to pretty print sigcleave results";

##
## More tests
##





