# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##
# $Id$

my $error;
use strict;
BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    $error = 0; 
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    use vars qw($TESTCOUNT);
    $TESTCOUNT = 38;
    plan tests => $TESTCOUNT;
}

use Bio::Tools::Est2Genome;
my $verbose = 0;

my $parser = new Bio::Tools::Est2Genome(-file   => Bio::Root::IO->catfile
					('t','data', 'hs_est.est2genome'));

ok($parser);
my $feature_set = $parser->parse_next_alignment;
ok(ref($feature_set), qr/ARRAY/i );

ok(scalar @$feature_set, 7);
my @exons = grep { $_->primary_tag eq 'Exon' } @$feature_set;
my @introns = grep { $_->primary_tag eq 'Intron' } @$feature_set;

my @expected_exons = ( [695,813,1,1,119,1],
		       [1377,1493,1,120,236,1],
		       [1789,1935,1,237,382,1],
		       [2084,2180,1,383,479,1]);
my @expected_introns = ( [814,1376,1],
			 [1494,1788,1],
			 [1936,2083,1] );

foreach my $e ( @exons ) {
    my $test_e = shift @expected_exons;
    my $i = 0;
    ok($e->query->start, $test_e->[$i++]);
    ok($e->query->end, $test_e->[$i++]);
    ok($e->query->strand, $test_e->[$i++]);

    ok($e->hit->start, $test_e->[$i++]);
    ok($e->hit->end, $test_e->[$i++]);
    ok($e->hit->strand, $test_e->[$i++]);
}
ok(! @expected_exons);
foreach my $intron ( @introns ) {
    my $test_i = shift @expected_introns;
    my $i = 0;
    ok($intron->start, $test_i->[$i++]);
    ok($intron->end, $test_i->[$i++]);
    ok($intron->strand, $test_i->[$i++]);
}
ok(! @expected_introns);
