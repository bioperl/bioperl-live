## Bioperl Test Harness Script for various modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..8\n"; 
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

$protein = 'VXCAAEFDFMEKETPLRYTKTLLLPVVLVVFVAIVRKIISDMWGVLAKQQTHVRKHQFDHGELVYHALQLLAYTALGILIMRLKLFLTPYMCVMASLICSRQLFGWLFCKVHPGAIVFVILAAMSIQGSANLQTQWKSTASLALET';

## need this to escape the strict() warnings...
$formatted_output = "";
my %results;


# Build object
test 2, $sigcleave_object = new Bio::Tools::Sigcleave(-ID   =>'test_sigcleave_seq',
                                                      -TYPE =>'amino',
			                              -SEQ  =>$protein); 
# Test raw result accessor
test 3, %results = $sigcleave_object->signals, "unable to get raw sigcleave results";

# Test formatted output method
test 4, $formatted_output = $sigcleave_object->pretty_print, "unable to pretty print sigcleave results";

##
## More tests
##

$sigcleave_object = new Bio::Tools::Sigcleave(-ID   =>'test_sigcleave_seq',
                                              -TYPE =>'amino',
		                              -SEQ  =>$protein); 
%results = $sigcleave_object->signals;

unless($results{111}) { print "not ok 5   expected to see a sigcleave score at position 111\n";
 } else { print "ok 5\n"; }

unless($results{130}) { print "not ok 6   expected to see a sigcleave score at position 130\n";
 } else { print "ok 6\n"; }

unless($results{111} == 3.7) { print "not ok 7   expected sigcleave score of 3.7 at position 111\n";
 } else { print "ok 7\n"; }

unless($results{130} == 4.9) { print "not ok 8   expected sigcleave score of 4.9 at position 130\n";
 } else { print "ok 8\n"; }





