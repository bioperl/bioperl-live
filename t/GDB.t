# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;    
    use vars qw($NUMTESTS);
    $NUMTESTS = 11;
    plan tests => $NUMTESTS;
    
    eval { require ('LWP/UserAgent.pm'); };
    if( $@ ) {
	print STDERR "Cannot load LWP::UserAgent, skipping tests\n";
	foreach ( 1..$NUMTESTS) { skip(1,1); }
	exit(0);
    }
}

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
ok $info->{'length'}, 226, 'value was '. $info->{'length'};

$marker = 'UT497';
$info = undef;
ok ($info = $gdb->get_info(-type=>'marker',
			     -id  => $marker));
ok $info->{gdbid}, 'GDB:198271', 'value was ' . $info->{gdbid};
ok $info->{primers}->[0], 'GGGTGACAGAACAAGACCT', 'value was ' . $info->{primers}->[0];
ok $info->{primers}->[1], 'ACCCATTAGCCTTGAACTGA', 'value was ' . $info->{primers}->[1];
ok $info->{'length'}, 155, 'value was '. $info->{'length'};
