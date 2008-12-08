# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 24);
	
	use_ok('Bio::Phenotype::MeSH::Term');
	use_ok('Bio::Phenotype::MeSH::Twig');
}
# For tests of Bio::DB::MeSH see t/DB.t

my $verbose = test_debug();

ok my $term = Bio::Phenotype::MeSH::Term->new(-verbose =>$verbose);
is $term->id('D000001'), 'D000001';
is $term->id, 'D000001';
is $term->name('Dietary Fats'), 'Dietary Fats';
is $term->name, 'Dietary Fats';
is $term->description('dietary fats are...'), 'dietary fats are...';
is $term->description, 'dietary fats are...';

ok my $twig = Bio::Phenotype::MeSH::Twig->new(-verbose =>$verbose);
is $twig->parent('Fats'), 'Fats';
is $twig->parent(), 'Fats';


ok $term->add_twig($twig);
is $term->each_twig(), 1;
is $twig->term, $term;

is $twig->add_sister('Bread', 'Candy', 'Cereals'), 3;
is $twig->add_sister('Condiments', 'Dairy Products'), 2;
is $twig->each_sister(), 5;
ok $twig->purge_sisters();
is $twig->each_sister(), 0;

is $twig->add_child('Butter', 'Margarine'), 2;
is $twig->each_child(), 2;
ok $twig->purge_children();
is $twig->each_child(), 0;
