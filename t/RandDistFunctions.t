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
    plan tests => 4;
}

use Bio::Tools::RandomDistFunctions;

my $dist = Bio::Tools::RandomDistFunctions->new;
ok(1);

ok($dist->rand_exponentional_distribution(1.0));
ok($dist->rand_geometric_distribution(0.9));
ok($dist->rand_normal_distribution);
