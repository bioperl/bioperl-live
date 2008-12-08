# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 3,
			   -requires_module => 'Bio::Ext::Align');
	
	use_ok('Bio::SearchDist');
}

my $dist = Bio::SearchDist->new();

isa_ok $dist, 'Bio::SearchDist';

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

is $dist->fit_Gaussian(1200), 1;

# TODO? there are no useful tests here!
