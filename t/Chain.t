# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
# Created: Wed Dec 13 15:52:33 GMT 2000
#
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..37\n";
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use lib '../';
use Bio::LiveSeq::Chain;

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

my $chain = Bio::LiveSeq::Chain::string2chain("abcdefghijklmnopqrstuvwxyz");
test 2, defined $chain && $chain ne '0';
test 3, Bio::LiveSeq::Chain::down_chain2string($chain) eq "abcdefghijklmnopqrstuvwxyz";
test 4, Bio::LiveSeq::Chain::down_chain2string($chain,undef,4) eq "abcd"; # default start=1
test 5, Bio::LiveSeq::Chain::down_chain2string($chain,1,4,6) eq "abcdef"; # last should override len
my $arrayref=Bio::LiveSeq::Chain::down_labels($chain,1,4);
test 6, $arrayref->[1] == 2;
$arrayref=Bio::LiveSeq::Chain::up_labels($chain,4,1);
test 7, $arrayref->[1] == 3;
$arrayref=Bio::LiveSeq::Chain::up_labels($chain);
test 8, scalar(@{$arrayref}) == 26; # total number of labels should be 26
test 9, Bio::LiveSeq::Chain::start($chain) eq '1';
test 10, Bio::LiveSeq::Chain::end($chain) eq '26';
test 11, Bio::LiveSeq::Chain::label_exists($chain,'4');
test 12, Bio::LiveSeq::Chain::label_exists($chain,'28') eq '0';
test 13, Bio::LiveSeq::Chain::down_get_pos_of_label($chain,4) eq '4';
test 14, Bio::LiveSeq::Chain::down_get_pos_of_label($chain,4,4) eq '1';
test 15, Bio::LiveSeq::Chain::up_get_pos_of_label($chain,26,1) eq '1';
test 16, Bio::LiveSeq::Chain::down_subchain_length($chain,1,4) eq '4';
test 17, Bio::LiveSeq::Chain::up_subchain_length($chain,4,1) eq '4';
test 18, Bio::LiveSeq::Chain::invert_chain($chain);
test 19, Bio::LiveSeq::Chain::invert_chain($chain);
test 20, Bio::LiveSeq::Chain::down_get_value_at_pos($chain,4) eq 'd';
test 21, Bio::LiveSeq::Chain::down_get_value_at_pos($chain,1,4) eq 'd';
test 22, Bio::LiveSeq::Chain::up_get_value_at_pos($chain,4) eq 'w';
test 23, Bio::LiveSeq::Chain::up_set_value_at_pos($chain,'W',4) &&  Bio::LiveSeq::Chain::up_get_value_at_pos($chain,4) eq 'W';
test 24, Bio::LiveSeq::Chain::down_set_value_at_pos($chain,'D',4) &&  Bio::LiveSeq::Chain::down_get_value_at_pos($chain,4) eq 'D';
test 25, Bio::LiveSeq::Chain::set_value_at_label($chain,'d',4) &&  Bio::LiveSeq::Chain::get_value_at_label($chain,4) eq 'd';
test 26, Bio::LiveSeq::Chain::down_get_label_at_pos($chain,1,4) eq '4';
test 27, Bio::LiveSeq::Chain::up_get_label_at_pos($chain,4) eq '23';
test 28, Bio::LiveSeq::Chain::is_downstream($chain,3,4);
test 29, Bio::LiveSeq::Chain::is_downstream($chain,4,3) eq '0';
test 30, Bio::LiveSeq::Chain::is_upstream($chain,4,3);
test 31, Bio::LiveSeq::Chain::is_upstream($chain,3,4) eq '0';
test 32, Bio::LiveSeq::Chain::splice_chain($chain,4,2) eq 'de';
test 33, Bio::LiveSeq::Chain::splice_chain($chain,7,undef,9) eq 'ghi';
my @array=Bio::LiveSeq::Chain::praeinsert_string($chain,"ghi",10);
test 34, $array[0] == 27 && $array[1] == 29;
my @array=Bio::LiveSeq::Chain::postinsert_string($chain,"de",3);
test 35, $array[0] == 30 && $array[1] == 31;
test 36, Bio::LiveSeq::Chain::up_chain2string($chain) eq "zyxWvutsrqponmlkjihgfedcba";
@array=Bio::LiveSeq::Chain::check_chain($chain);
test 37, $array[0] == 1 && $array[1] == 1 && $array[2] == 1 && $array[3] == 1;

