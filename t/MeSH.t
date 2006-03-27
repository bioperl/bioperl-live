# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

my $error;

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
	$NUMTESTS = 23;

	plan tests => $NUMTESTS;

	eval { require IO::String; 
			 require LWP::UserAgent;
			 require HTTP::Request::Common;
       };
	if( $@ ) {
		print STDERR "IO::String or LWP::UserAgent or HTTP::Request not installed. This means the MeSH modules are not usable. Skipping tests.\n";
		for( 1..$NUMTESTS ) {
			skip("IO::String, LWP::UserAgent,or HTTP::Request not installed",1);
		}
		$error = 1;
	}
}
# For tests of Bio::DB::MeSH see t/DB.t

if( $error ==  1 ) {
    exit(0);
}

END {
	foreach ( $Test::ntest..$NUMTESTS) {
		skip('unable to run all of the MeSH.t tests, skipping',1);
    }
}

require Bio::Phenotype::MeSH::Term;
require Bio::Phenotype::MeSH::Twig;
require Bio::DB::MeSH;
ok 1;

my $verbose = 0;

ok my $term = Bio::Phenotype::MeSH::Term->new(-verbose =>$verbose);
ok $term->id('D000001'), 'D000001';
ok $term->id, 'D000001';
ok $term->name('Dietary Fats'), 'Dietary Fats';
ok $term->name, 'Dietary Fats';
ok $term->description('dietary fats are...'), 'dietary fats are...';
ok $term->description, 'dietary fats are...';

ok my $twig = Bio::Phenotype::MeSH::Twig->new(-verbose =>$verbose);
ok $twig->parent('Fats'), 'Fats';
ok $twig->parent(), 'Fats';


ok $term->add_twig($twig);
ok $term->each_twig(), 1;
ok $twig->term, $term;

ok $twig->add_sister('Bread', 'Candy', 'Cereals'), 3;
ok $twig->add_sister('Condiments', 'Dairy Products'), 2;
ok $twig->each_sister(), 5;
ok $twig->purge_sisters();
ok $twig->each_sister(), 0;

ok $twig->add_child('Butter', 'Margarine'), 2;
ok $twig->each_child(), 2;
ok $twig->purge_children();
ok $twig->each_child(), 0;


