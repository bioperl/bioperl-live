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

    $NUMTESTS = 26;
    plan tests => $NUMTESTS;

}

END {
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('unable to run all of the tests depending on web access',1);
    }
}


#use Data::Dumper;
use Bio::Phenotype::MeSH::Term;
use Bio::Phenotype::MeSH::Twig;
use Bio::DB::MeSH;
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


eval {
    ok my $mesh = new Bio::DB::MeSH(-verbose => $verbose);
    ok my $t=$mesh->get_exact_term('Butter');
    ok $t->each_twig(), 3;
    #print Dumper $t;
};
