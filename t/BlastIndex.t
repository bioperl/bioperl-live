# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

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
    plan tests => 13;
}

use Bio::Tools::BPlite;
use Bio::Index::Blast;
use Bio::Root::IO;

END {  unlink qw( Wibbl ); }

ok(1);

my $index = new Bio::Index::Blast(-filename => 'Wibbl',
				  -write_flag => 1);
ok($index);
$index->make_index(Bio::Root::IO->catfile("t","data","multi_blast.bls"));
ok(-e "Wibbl");

foreach my $id ( qw(CATH_RAT PAPA_CARPA) ) {

    my $fh = $index->get_stream($id);
    ok($fh);
    ok( ! eof($fh) );
    my $report = new Bio::Tools::BPlite(-fh => $fh);
    ok($report->query, qr/$id/);
    ok( $report->nextSbjct);
    ok( $index->fetch_report($id)->query, qr/$id/);
}

