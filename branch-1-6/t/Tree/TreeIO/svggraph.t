# -*-Perl-*- Test Harness script for Bioperl
# $Id: TreeIO.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 4);
    use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

my $treeio = Bio::TreeIO->new(-verbose => $verbose,
			  -file   => test_input_file('test.nhx'));
my $tree = $treeio->next_tree;

SKIP: {
	test_skip(-tests => 3, -requires_module => 'SVG::Graph');
	my $FILE3 = test_output_file();
	my $treeout3 = Bio::TreeIO->new(-format => 'svggraph',
					-file => ">$FILE3");
	ok($treeout3);
	eval {$treeout3->write_tree($tree);};
	ok(!$@);
	ok (-s $FILE3);
}
