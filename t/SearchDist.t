# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

use Test;
use strict;
BEGIN { 
    eval { local * STDERR; require Bio::Ext::Align };
    if ( $@ ) {
	plan tests => 1;
	skip(1, 1, 1,'Bio::Ext::Align not loaded');
	print STDERR "\tBio::Ext::Align not loaded\n";
	exit(0);
    }
    plan tests => 3;
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

    
ok $dist->fit_evd(), 1;

foreach my $score ( @scores ) {
    my $evalue = $dist->evalue($score);
    print SDTERR "Score $score had an evalue of $evalue\n";
}
