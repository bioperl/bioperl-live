# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

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
    plan tests => 3;

    eval { local * STDERR; require Bio::Ext::Align };
    if ( $@ ) {
	foreach ( 1..3) {
	    skip(1,'Bio::Ext::Align not loaded');
	}
	exit(0);
    }
}
use Bio::SearchDist;
ok(1);

my $dist = new Bio::SearchDist;

ok ref($dist), 'Bio::SearchDist';

my @scores;
foreach my $i ( 1..5000 ) {
    my $score = rand(1300);
    #print STDERR "Got $score\n";
    $dist->add_score($score);
    push(@scores,$score);
}

# this just checks that this routine runs ;)
# as the distribution is not gaussian, it gives
# non-sensical results    

ok $dist->fit_Gaussian(1200), 1;

