# -*-Perl-*-
# $Id$
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
    plan tests => 11;
}

use Bio::PopGen::Simulation::Coalescent;
use Bio::PopGen::Statistics;
use Bio::TreeIO;
ok(1);

use vars qw($FILE1);
$FILE1 = 'out.tre';
END { unlink $FILE1; }
 
my $ssize = 5;
my $sim = new Bio::PopGen::Simulation::Coalescent(-sample_size => $ssize);
my $stats = Bio::PopGen::Statistics->new();
my $tree = $sim->next_tree;

ok($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = new Bio::TreeIO(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
my $mutcount = 100;
$sim->add_Mutations($tree, $mutcount);

my $leaves = [$tree->get_leaf_nodes];
# $stats->verbose(1);
my $pi = $stats->pi($leaves);
ok( $pi > 0 , 1, 'pi');


# theta is num seg sites / sum(1/(numsams-1))
my $theta = $stats->theta(scalar @$leaves, $mutcount);
ok($theta,48, 'theta');

my $tD = Bio::PopGen::Statistics->tajima_D($leaves);
ok(defined $tD,1, 'tajimaD');

my $seg_sites = Bio::PopGen::Statistics->segregating_sites_count($leaves);
ok($seg_sites,$mutcount,
   'all the mutations should be polymorphic (by definition)');
my $single = Bio::PopGen::Statistics->singleton_count($leaves);
my $flD = Bio::PopGen::Statistics->fu_and_li_D($leaves,$single);
ok(defined $flD,1,'fu and li D');

my $flD_star = $stats->fu_and_li_D_star($leaves);
ok(defined $flD_star,1,'fu and li D*');

my $flF = $stats->fu_and_li_F($leaves,$single);
ok(defined $flF, 1,'fu and li F');

my $flFstar = $stats->fu_and_li_F_star($leaves);
ok(defined $flF, 1,'fu and li F');
