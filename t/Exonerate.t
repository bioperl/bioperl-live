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
    $NTESTS = 45;
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
my @data = ( [qw(ln27 416 Contig124 939
		 293 416 -1 
		 1   124 1 
		 
		 107 292 -1 
		 178 363 1 
		 
		 66 106 -1
		 899 939 1
		 )],
	     [qw(ln74 644 Contig275 1296 
		 601 644 -1
		 901 944 1
		 
		 436 600 -1
		 998 1162    1

		 386 435 -1
		 1247 1296 1
		 
		 )] );

while( my $r = $searchio->next_result ) {
    my $d = shift @data;
    ok($r->query_name, shift @$d);
    skip( 'no query length available in default output',
	  $r->query_length, shift @$d);
    my $h = $r->next_hit;
    ok($h->name, shift @$d);
    skip( 'no hit length available in default output',$h->length, shift @$d);
    while( my $hsp = $h->next_hsp ) {
	ok($hsp->query->start, shift @$d);
	ok($hsp->query->end, shift @$d);
	ok($hsp->query->strand, shift @$d);
	
	ok($hsp->hit->start, shift @$d);
	ok($hsp->hit->end, shift @$d);
	ok($hsp->hit->strand, shift @$d);
    }
}
