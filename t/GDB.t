# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##


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
use Test;
use strict;

BEGIN { plan tests => 11 }


use Bio::DB::GDB;
my $verbose = 0;

my ($gdb, $marker, $info);
# get a single seq
ok defined ( $gdb = new Bio::DB::GDB(-verbose=>$verbose) );     
$marker = 'D1S234';
ok ($info = $gdb->get_info(-type=>'marker',
			   -id  => $marker));
ok $info->{gdbid}, 'GDB:188296', 'value was ' . $info->{gdbid};
ok $info->{primers}->[0], 'GCCCAGGAGGTTGAGG', 'value was ' . $info->{primers}->[0];
ok $info->{primers}->[1], 'AAGGCAGGCTTGAATTACAG', 'value was ' . $info->{primers}->[1];
ok $info->{length}, 226, 'value was '. $info->{length};

$marker = 'UT497';
$info = undef;
ok ($info = $gdb->get_info(-type=>'marker',
			     -id  => $marker));
ok $info->{gdbid}, 'GDB:198271', 'value was ' . $info->{gdbid};
ok $info->{primers}->[0], 'GGGTGACAGAACAAGACCT', 'value was ' . $info->{primers}->[0];
ok $info->{primers}->[1], 'ACCCATTAGCCTTGAACTGA', 'value was ' . $info->{primers}->[1];
ok $info->{length}, 155, 'value was '. $info->{length};
