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
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use vars qw($NTESTS);
    $NTESTS = 49;
    $error = 0;
    use Test::More;
    plan tests => $NTESTS;
	use_ok('Bio::SearchIO');
	use_ok('Bio::Root::IO');
	use_ok('Bio::SearchIO::Writer::HitTableWriter');
	use_ok('Bio::SearchIO::Writer::HTMLResultWriter');
}

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


my $val;

while( my $r = $searchio->next_result ) {
    my $d = shift @data;
    is($r->query_name, shift @$d);
	SKIP: {
		$val = shift @$d;
		skip('no query length available in default output',1);
	    is($r->query_length, $val);
		   };
    
    my $h = $r->next_hit;
    is($h->name, shift @$d);
	SKIP: {
		$val = shift @$d;
		skip( 'no hit length available in default output',1);
	    is($h->length, $val);
		   };
    while( my $hsp = $h->next_hsp ) {
	is($hsp->query->start, shift @$d);
	is($hsp->query->end, shift @$d);
	is($hsp->query->strand, shift @$d);
	
	is($hsp->hit->start, shift @$d);
	is($hsp->hit->end, shift @$d);
	is($hsp->hit->strand, shift @$d);
    }
}
