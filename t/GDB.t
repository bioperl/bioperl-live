# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
my $error;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    $error = 0;
    use Test;    
    use vars qw($NUMTESTS);
    $NUMTESTS = 11;
    plan tests => $NUMTESTS;
    
    eval { require ('LWP/UserAgent.pm'); require('HTML/Parser.pm');	
	   require ('HTTP/Request/Common.pm');
	 };
    if( $@ ) {
	print STDERR "Cannot load LWP::UserAgent or HTML::Parser, skipping tests\n";
	foreach ( 1..$NUMTESTS) { skip(1,1); }
	$error = 1;
    } 
    if( $] < 5.005 ) {
	print STDERR "GDB parsing does not work with 5.005 or lower Perl versions.\n";
	foreach ( 1..$NUMTESTS) { skip(1,1); }
	$error = 1;
    }
}

if( $error == 1 ) {
    exit(0);
}

require Bio::DB::GDB;
my $verbose = 0;

my ($gdb, $marker, $info);
# get a single seq
$marker = 'D1S234';
eval { 
    ok defined ( $gdb = new Bio::DB::GDB(-verbose=>$verbose) );     
    ok($info = $gdb->get_info(-type=>'marker',
			      -id  => $marker));
};
if( $@ ) {
    warn "Warning: Couldn't connect to GDB website with Bio::DB::GDB.pm!\nError: Do you have network access? Skipping all other tests";
    foreach ( $Test::ntest..$NUMTESTS ) { skip(1,1, 'no network access'); }
    exit(0);
}


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
