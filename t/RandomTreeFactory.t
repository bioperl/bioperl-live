# -*-Perl-*-
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
    plan tests => 4;
}

use Bio::Tree::RandomFactory;
use Bio::TreeIO;
use Bio::Tree::Statistics;
ok(1);

use vars qw($FILE1);
$FILE1 = 'out.tre';
END { unlink $FILE1; }
 
my $ssize = 10;
my $factory = new Bio::Tree::RandomFactory(-sample_size => $ssize);
my $stats = new Bio::Tree::Statistics();

my $tree = $factory->next_tree;

ok($tree->get_nodes, ($ssize * 2 - 1));

my $treeio = new Bio::TreeIO(-format => 'newick', -file => ">$FILE1");

$treeio->write_tree($tree);
undef $treeio;

ok(-s $FILE1);
my $mutcount = 100;
$factory->add_Mutations($tree, $mutcount);

my $flD = $stats->fu_and_li_D($tree);

ok(defined $flD,1, 'fu and li D');

my $flD_star = $stats->fu_and_li_D_star($tree);
ok(defined $flD_star,1, 'fu and li D*');

my $flF = $stats->fu_and_li_F($tree);
ok(defined $flF,1,'fu and li F');

my $tD = $stats->tajima_D($tree);
ok(defined $tD,1, 'tajimaD');

my $theta = $stats->theta($tree);
ok(defined $theta,1, 'theta');

$stats->verbose(1);
my $pi = $stats->pi($tree);
#ok($pi,1, 'pi');
