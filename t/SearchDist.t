# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use Test;
use strict;
BEGIN { 
    plan tests => 5; 
    eval {
	require Bio::SearchDist;
	require Bio::Ext::Align;
    };
    if ( $@ ) {
	for(my $i=1;$i<=5;$i++) {	    
	    skip(1, 'could not find Bio::SearchDist');
	}
	exit(0);
    }
}

ok(1);

my $dist = new Bio::SearchDist;

ok ref($dist), 'Bio::SearchDist';
my @scores =  qw( 100 200 120 121 78 165 215 6 18);
foreach my $score ( @scores ) {
    $dist->add_score($score);
}
    
ok $dis->fit_evd(), 1;

foreach $score ( @scores ) {
    my $evalue = $dist->evalue($score);
    print "Score $score had an evalue of $evalue\n";
}
