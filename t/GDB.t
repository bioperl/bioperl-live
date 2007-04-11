# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    $NUMTESTS = 12;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't/lib';
	}
	use Test::More;
	
	eval {
		require LWP::UserAgent;
        require HTML::Parser;
        require HTTP::Request::Common;
	};
	if ($@) {
		plan skip_all => 'HTML::Parser or LWP::UserAgent or HTTP::Request not installed. This means Bio::DB::GDB modules is not usable. Skipping tests';
	}
    elsif (!$DEBUG) {
		plan skip_all => 'Skipping all tests since they require network access, set BIOPERLDEBUG=1 to test';
	}
	else {
		plan tests => $NUMTESTS;
	}
    use_ok('Bio::DB::GDB');
}

my $verbose = -1;

my ($gdb, $marker, $info);
# get a single seq

SKIP: {
    $marker = 'D1S234';

    $gdb = new Bio::DB::GDB(-verbose=>$verbose);
    
    eval {
        $info = $gdb->get_info(-type=>'marker',
                      -id  => $marker);
    };
    
    if( $@ || ! defined $info) {
        skip("Warning: Couldn't connect to GDB website!  Skipping all other tests",11);
    }
    
    ok $gdb;
    ok $info;
    is $info->{gdbid}, 'GDB:188296', 'info was ' . $info->{gdbid};
    is $info->{primers}->[0], 'GCCCAGGAGGTTGAGG', 'info was ' . $info->{primers}->[0];
    is $info->{primers}->[1], 'AAGGCAGGCTTGAATTACAG', 'info was ' . $info->{primers}->[1];
    is $info->{'length'}, 226, 'info was '. $info->{'length'};
    
    $marker = 'UT497';
    $info = undef;
    eval { 
        $info = $gdb->get_info(-type=>'marker',
                         -id  => $marker);
    };
    if( $@ || ! defined $info) {
        skip("Warning: Couldn't connect to GDB website! Skipping all other tests",5);
    }
    ok $info;
    is $info->{gdbid}, 'GDB:198271', 'info was ' . $info->{gdbid};
    is $info->{primers}->[0], 'GGGTGACAGAACAAGACCT', 'info was ' . $info->{primers}->[0];
    is $info->{primers}->[1], 'ACCCATTAGCCTTGAACTGA', 'info was ' . $info->{primers}->[1];
    is $info->{'length'}, 155, 'info was '. $info->{'length'};
}

