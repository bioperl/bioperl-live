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
    plan tests => 3;
}

use Bio::TreeIO;
my $verbose = 0;

my $treeio = new Bio::TreeIO(-verbose => $verbose,
			     -format => 'nhx',
			     -file   => Bio::Root::IO->catfile('t','data', 
							       'test.nhx'));

my $tree = $treeio->next_tree;

my @nodes = $tree->find_node('ADH2');
ok(@nodes, 2);

if( $verbose ) {
    print "nodes are: \n",
    join(", ", map {  $_->id . ":". (defined $_->branch_length ? 
				     $_->branch_length : '' ) } @nodes), "\n";
}

$treeio = new Bio::TreeIO(-format => 'newick',
			  -file   => Bio::Root::IO->catfile('t','data',
							    'test.nh'));
$tree = $treeio->next_tree;
if( $verbose ) { 
    my $out = new Bio::TreeIO(-format => 'tabtree');
    
    $out->write_tree($tree);
}

my @hADH = ( $tree->find_node('hADH1'),
	     $tree->find_node('hADH2') );
my ($n4) = $tree->find_node('yADH4');

ok($tree->is_monophyletic(-nodes    => \@hADH,
			  -outgroup => $n4));

my @mixgroup = ( $tree->find_node('hADH1'),
		 $tree->find_node('yADH2'),
		 $tree->find_node('yADH3'),
		 );

my ($iADHX) = $tree->find_node('iADHX');

ok(! $tree->is_monophyletic(-nodes   => \@mixgroup,
			    -outgroup=> $iADHX));
