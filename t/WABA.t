# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error;

use strict;
use lib '.';

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use vars qw($NTESTS);
    $NTESTS = 62;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 
}

if( $error == 1 ) {
    exit(0);
}

use Bio::SearchIO;
use Bio::Root::IO;

my $wabain = new Bio::SearchIO(-format => 'waba',
			       -file   => Bio::Root::IO->catfile('t','data',
								 'test.waba'));

# These won't look the same as the WABA file because Jim's code is 0 based
# while we (bioperl) are 1 based.
my @results = ( 
		[ 'U57623', 'pair1_hs.fa', 'pair1_mm.fa',
		  [ 'U02884', 3, 
		    [qw(3833 34 2972 1 243 3688 1 40.9)],
		    [qw(4211 3022 6914 1 3705 6848 1 43.7)],
		    [qw(2218 7004 9171 1 6892 8712 1 50.3)],
		    ], 
		  ],
		[ 'X57152', 'pair9_hs.fa', 'pair9_mm.fa',
		  [ 'X80685', 1, 
		    [qw(7572 4 5845 1 632 7368 1 46.8)],
		    ], 
		  ]
		);
while( my $wabar = $wabain->next_result )  {
    my @r = @{shift @results};
    ok($wabar->query_name, shift @r);
    ok($wabar->query_database, shift @r);
    ok($wabar->database_name, shift @r);
    while( my $wabah = $wabar->next_hit ) {
	my (@h) = @{shift @r};
	ok( $wabah->name, shift @h);
	ok( $wabah->hsps(), shift @h);
	while( my $wabahsp = $wabah->next_hsp  ) {
	    my ( @hsp) = @{shift @h};
	    ok($wabahsp->length('total'), shift @hsp);
	    ok($wabahsp->query->start, shift @hsp);
	    ok($wabahsp->query->end, shift @hsp);
	    ok($wabahsp->strand('query'), shift @hsp);
	    ok($wabahsp->start('hit'), shift @hsp);
	    ok($wabahsp->end('subject'), shift @hsp);
	    ok($wabahsp->subject->strand, shift @hsp);
	    ok(length($wabahsp->query_string), $wabahsp->length('total'));
	    ok(length($wabahsp->hit_string), $wabahsp->length('total'));
	    ok(length($wabahsp->hmmstate_string), $wabahsp->length('total'));
	    my $hs = $wabahsp->hit_string;
	    ok($wabahsp->gaps('hit'), $hs  =~ tr/\-//);
	    my $qs = $wabahsp->query_string;
	    ok($wabahsp->gaps('query'),  $qs =~ tr/\-//);
	    ok(sprintf("%.1f",$wabahsp->percent_identity),shift @hsp);
	}
    }
}
