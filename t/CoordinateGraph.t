#-*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 7;
	use_ok('Bio::Coordinate::Graph');
}

ok my $graph = Bio::Coordinate::Graph->new();

# graph structure
my $dag = {
	   9  => [],
	   8  => [9],
	   7  => [],
	   6  => [7, 8],
	   5  => [],
	   4  => [5],
	   3  => [6],
	   2  => [3, 4, 6],
	   1  => [2]
	  };

ok $graph->hash_of_arrays($dag);


my $a = 1;
my $b = 6;
is my @a = $graph->shortest_path($a, $b), 3;
#print join (", ", @a), "\n";

$a = 7;
$b = 8;
is @a = $graph->shortest_path($a, $b), 1;


$a = 8;
$b = 9;
is @a = $graph->shortest_path($a, $b), 2;
$b = 2;
is @a = $graph->shortest_path($a, $b), 3;


