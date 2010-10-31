# -*-Perl-*- Test Harness script for Bioperl

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 19);

    use_ok('Bio::SearchIO');
}

my $axtparser = Bio::SearchIO->new(-format => 'axt',
				  -file   => test_input_file('test_data.axt'));

my $result = $axtparser->next_result;
is($result->query_name, 'chr19');

my $hit = $result->next_hit;
is($hit->name, 'chr11');

my $hsp    = $hit->next_hsp;
is($hsp->query->start,3001012);
is($hsp->query->end,3001075);
is($hsp->query->length, 64);
is($hsp->query->strand, 1);

is($hsp->hit->start,70568380);
is($hsp->hit->end,70568443);
is($hsp->hit->length, 64);
is($hsp->hit->strand, -1);

# next HSP returns next alignment, but shouldn't this be the next_hit instead????
# what happens if the hit name changes?

$hsp    = $hit->next_hsp;
is($hsp->query->start,3008279);
is($hsp->query->end,3008357);
is($hsp->query->length, 79);
is($hsp->query->strand, 1);

is($hsp->hit->start,70573976);
is($hsp->hit->end,70574054);
is($hsp->hit->length, 79);
is($hsp->hit->strand, -1);

