# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $error);


BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 17;
    plan tests => $NUMTESTS;
    eval { require IO::String; 
	   require Bio::Tools::Phylo::Molphy;}; 
    if( $@ ) { print STDERR "no IO string installed\n"; 
	       $error = 1;
	   }
}

END { 
    foreach ( $Test::ntest .. $NUMTESTS ) {
	skip("unable to run all of the Molphy tests",1);
    }
}


exit(0) if( $error );

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $inmolphy = new Bio::Tools::Phylo::Molphy(-file => 't/data/lysozyme6.simple.protml');
ok($inmolphy);
my $r = $inmolphy->next_result;
ok($r);
ok($r->model, 'JTT');
ok($r->search_space,50);
my @trees;
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
ok(@trees,5);
 $inmolphy = new Bio::Tools::Phylo::Molphy(-file => 't/data/lysozyme6.protml');
ok($inmolphy);
$r = $inmolphy->next_result;
ok($r->model, 'JTT');
ok($r->search_space,50);
@trees = ();
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
ok(@trees,5);

ok($trees[0]->score, -1047.8);
ok($trees[-1]->id, 9);

my $tpm = $r->transition_probability_matrix;
ok($tpm->{'Val'}->{'Val'}, -122884);
ok($tpm->{'Ala'}->{'Arg'}, 2710);

my $sub_mat = $r->substitution_matrix;
ok($sub_mat->{'Val'}->{'Tyr'}, 50);
ok($sub_mat->{'Arg'}->{'Ile'}, 72);
ok($sub_mat->{'Met'}->{'Met'}, '');

my %fmat = $r->residue_frequencies();
ok($fmat{'D'}->[0], 0.052);
