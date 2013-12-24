use strict;
use warnings;

BEGIN {
    use Bio::Root::Test;

    test_begin(-tests => 7);

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

$a = 7;
$b = 8;
is @a = $graph->shortest_path($a, $b), 1;

$a = 8;
$b = 9;
is @a = $graph->shortest_path($a, $b), 2;
$b = 2;
is @a = $graph->shortest_path($a, $b), 3;
