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
    $NTESTS = 17;
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
my @data = ( [qw(Contig124 ln27 575 764 -1 835 1024 1)],
	     [qw(Contig124 ln27 810 940 -1 1017 1147 1)] );
while( my $r = $searchio->next_result ) {
    my $d = shift @data;
    ok($r->query_name, shift @$d);
    my $h = $r->next_hit;
    ok($h->name, shift @$d);
    my $hsp = $h->next_hsp;
    ok($hsp->query->start, shift @$d);
    ok($hsp->query->end, shift @$d);
    ok($hsp->query->strand, shift @$d);

    ok($hsp->hit->start, shift @$d);
    ok($hsp->hit->end, shift @$d);
    ok($hsp->hit->strand, shift @$d);
}
