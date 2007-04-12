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
	eval { require Test::More; };
	$error = 0;
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;

	$NUMTESTS = 18;
	eval { require IO::String;
			};
	if( $@ ) {
		plan skip_all => 'IO::String not installed. Skipping Molphy tests'
	} else {
		plan tests => $NUMTESTS;
	}
	use_ok 'Bio::Tools::Phylo::Molphy';
}

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
is($r->model, 'JTT');
is($r->search_space,50);
my @trees;
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
is(@trees,5);
 $inmolphy = new Bio::Tools::Phylo::Molphy(-file => 't/data/lysozyme6.protml');
ok($inmolphy);
$r = $inmolphy->next_result;
is($r->model, 'JTT');
is($r->search_space,50);
@trees = ();
while( my $t = $r->next_tree ) { 
    push @trees, $t;
}
is(@trees,5);

is($trees[0]->score, -1047.8);
is($trees[-1]->id, 9);

my $tpm = $r->transition_probability_matrix;
is($tpm->{'Val'}->{'Val'}, -122884);
is($tpm->{'Ala'}->{'Arg'}, 2710);

my $sub_mat = $r->substitution_matrix;
is($sub_mat->{'Val'}->{'Tyr'}, 50);
is($sub_mat->{'Arg'}->{'Ile'}, 72);
is($sub_mat->{'Met'}->{'Met'}, '');

my %fmat = $r->residue_frequencies();
is($fmat{'D'}->[0], 0.052);
