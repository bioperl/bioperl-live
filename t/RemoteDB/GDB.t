# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 12,
			   -requires_modules => [qw(LWP::UserAgent
									    HTML::Parser
										HTTP::Request::Common)],
			   -requires_networking => 1);
	
    use_ok('Bio::DB::GDB');
}

my $verbose = test_debug();

my ($gdb, $marker, $info);
# get a single seq

SKIP: {
    $marker = 'D1S234';

    $gdb = Bio::DB::GDB->new(-verbose=>$verbose);
    
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
