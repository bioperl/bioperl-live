# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

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
    plan tests => 3;

    eval { local * STDERR; require Bio::Ext::Align };
    if ( $@ ) {
	foreach ( 1..3) {
	    skip('Bio::Ext::Align not loaded',1);
	}
        $error = 1;
    }
}

if( $error == 1 ) { 
    exit(0);
}

require Bio::SearchDist;
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

