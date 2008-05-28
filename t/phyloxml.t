# -*-Perl-*- Test Harness script for Bioperl
# $Id: phyloxml.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
	use lib 't/lib';
  use BioperlTest;
    
  test_begin(-tests => 4);
	
	use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

ok my $treeio = Bio::TreeIO->new(
           -verbose => $verbose,
			     -format => 'phyloxml',
			     -file   => test_input_file('phyloxml_small.xml'));

my $tree = $treeio->next_tree;
isa_ok($tree, 'Bio::Tree::TreeI');

TODO: {
  local $TODO = 'write_tree not implemented yet';
  my $FILE1 = test_output_file();
  $treeio = Bio::TreeIO->new(-verbose => $verbose,
          -format => 'phyloxml',
          -file   => ">$FILE1");
  $treeio->write_tree($tree);
  undef $treeio;
  ok( -s $FILE1 );
}

