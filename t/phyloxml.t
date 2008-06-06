# -*-Perl-*- Test Harness script for Bioperl
# $Id: phyloxml.t 14580 2008-03-01 17:01:30Z cjfields $

use strict;

BEGIN {
	use lib 't/lib';
  use BioperlTest;
  use XML::LibXML;
  if (1000*$] < 5008) {
     plan skip_all => "Reader interface only supported in Perl >= 5.8";
     exit;
  } elsif (XML::LibXML::LIBXML_VERSION() <= 20620) {
     plan skip_all => "Reader not supported for libxml2 <= 2.6.20";
     exit;
  } else {
    test_begin(-tests => 18);
  } 
	
	use_ok('Bio::TreeIO');
}

my $verbose = test_debug();

ok my $treeio = Bio::TreeIO->new(
           -verbose => $verbose,
			     -format => 'phyloxml',
			     -file   => test_input_file('phyloxml_small.xml'));

my $tree = $treeio->next_tree;
isa_ok($tree, 'Bio::Tree::TreeI');

ok my $treeio = Bio::TreeIO->new(
           -verbose => $verbose,
			     -format => 'phyloxml',
			     -file   => test_input_file('phyloxml_examples.xml'));

while ( $tree = $treeio->next_tree ) {
  isa_ok($tree, 'Bio::Tree::TreeI');
}

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

