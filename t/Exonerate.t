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
    $NTESTS = 41;
    $error = 0;

    use Test;
    plan tests => $NTESTS; 
}

if( $error == 1 ) {
    exit(0);
}

use Bio::SearchIO;
use Bio::Root::IO;
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SearchIO::Writer::HTMLResultWriter;

ok(1);
my ($searchio, $result,$hit,$hsp);

$searchio = new Bio::SearchIO(-file => 
			      Bio::Root::IO->catfile(qw(t data 
							testdat.exonerate)),
			      -format => 'exonerate');
my @data = ( [qw(ln27 Contig124 
		 292 417 -1 
		 1 125 1 
		 
		 106 293 -1 
		 178 364 1 
		 
		 66 107 -1
		 899 940 1
		 )],
	     [qw(ln74 Contig275 
		 600 645 -1
		 901 945 1
		 
		 435 601 -1
		 998 1163 1
		 
		 386 436 -1
		 1247 1297 1
		 )] );

while( my $r = $searchio->next_result ) {
    my $d = shift @data;
    ok($r->query_name, shift @$d);
    my $h = $r->next_hit;
    ok($h->name, shift @$d);
    while( my $hsp = $h->next_hsp ) {
	ok($hsp->query->start, shift @$d);
	ok($hsp->query->end, shift @$d);
	ok($hsp->query->strand, shift @$d);
	
	ok($hsp->hit->start, shift @$d);
	ok($hsp->hit->end, shift @$d);
	ok($hsp->hit->strand, shift @$d);
    }
}
