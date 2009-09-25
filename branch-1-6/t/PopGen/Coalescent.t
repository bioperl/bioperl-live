# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 13);
	
	use_ok('Bio::PopGen::Simulation::Coalescent');
	use_ok('Bio::PopGen::Statistics');
	use_ok('Bio::TreeIO');
}

my $ssize = 5;
my $sim = Bio::PopGen::Simulation::Coalescent->new(-sample_size => $ssize);
my $stats = Bio::PopGen::Statistics->new();
my $tree = $sim->next_tree;

is($tree->get_nodes, ($ssize * 2 - 1));

my $FILE1 = test_output_file();
my $treeio = Bio::TreeIO->new(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
my $mutcount = 100;
$sim->add_Mutations($tree, $mutcount);

my $leaves = [$tree->get_leaf_nodes];
# $stats->verbose(1);
my $pi = $stats->pi($leaves);
is( $pi > 0 , 1, 'pi');


# theta is num seg sites / sum(1/(numsams-1))
my $theta = $stats->theta(scalar @$leaves, $mutcount);
is($theta,48, 'theta');

my $tD = Bio::PopGen::Statistics->tajima_D($leaves);
is(defined $tD,1, 'tajimaD');

my $seg_sites = Bio::PopGen::Statistics->segregating_sites_count($leaves);
is($seg_sites,$mutcount,
   'all the mutations should be polymorphic (by definition)');
my $single = Bio::PopGen::Statistics->singleton_count($leaves);
my $flD = Bio::PopGen::Statistics->fu_and_li_D($leaves,$single);
is(defined $flD,1,'fu and li D');

my $flD_star = $stats->fu_and_li_D_star($leaves);
is(defined $flD_star,1,'fu and li D*');

my $flF = $stats->fu_and_li_F($leaves,$single);
is(defined $flF, 1,'fu and li F');

my $flFstar = $stats->fu_and_li_F_star($leaves);
is(defined $flF, 1,'fu and li F');
