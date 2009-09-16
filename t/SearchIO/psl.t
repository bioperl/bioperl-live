# -*-Perl-*- Test Harness script for Bioperl
# $Id: UCSCParsers.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 53);

    use_ok('Bio::SearchIO');
}

my $pslparser = Bio::SearchIO->new(-format => 'psl',
				  -file   => test_input_file('sbay_c545-yeast.BLASTZ.PSL'));

my $result = $pslparser->next_result;
is($result->query_name, 'I');
is($result->query_length, 230203);

my $hit    = $result->next_hit;
is($hit->name, 'sbay_c545');
is($hit->length, 28791);
my $hsp    = $hit->next_hsp;
is($hsp->query->start,139871);
is($hsp->query->end,141472);
is($hsp->query->length, 1602);
is($hsp->query->strand, 1);
is($hsp->hit->strand, 1);
my $q_gapblocks = $hsp->gap_blocks('query');
is(scalar @$q_gapblocks, 24);
is($q_gapblocks->[0]->[1],45);
is($q_gapblocks->[1]->[1],10);
is($q_gapblocks->[1]->[0],139921);

$hsp       = $hit->next_hsp;
$hsp       = $hit->next_hsp;
is($hsp->hit->start,27302);
is($hsp->hit->end,27468);
is($hsp->hit->length,167);
is($hsp->query->start, 123814);
is($hsp->query->end, 123972);
is($hsp->query->length, 159);
is($hsp->query->strand,-1);

$q_gapblocks = $hsp->gap_blocks('query');
is(scalar @$q_gapblocks, 4);
is($q_gapblocks->[0]->[1],116);
is($q_gapblocks->[1]->[1],4);
is($q_gapblocks->[1]->[0],123856);

#-----------------------------------

$pslparser = Bio::SearchIO->new(-format => 'psl',
			       -file   => test_input_file('blat.psLayout3'));

$result = $pslparser->next_result;
is($result->query_name, 'sequence_10');
is($result->query_length, 1775);

$hit    = $result->next_hit;
is($hit->name, 'sequence_10');
is($hit->length, 1775);
$hsp    = $hit->next_hsp;
is($hsp->query->start,1);
is($hsp->query->end,1775);
is($hsp->query->length,1775);
is($hsp->query->strand,1);
is($hsp->hit->strand,1);
$q_gapblocks = $hsp->gap_blocks('query');
is(scalar @$q_gapblocks, 1);
is($q_gapblocks->[0]->[1],1775);
is($q_gapblocks->[1]->[1],undef);
is($q_gapblocks->[1]->[0],undef);

$hsp       = $hit->next_hsp;
is($hsp->hit->start,841);
is($hsp->hit->end,1244);
is($hsp->query->start, 841);
is($hsp->query->end, 1244);
is($hsp->query->length, 404);
is($hsp->query->strand,-1);
is($hsp->hit->strand,1);

$q_gapblocks = $hsp->gap_blocks('query');
is(scalar @$q_gapblocks, 4);
is($q_gapblocks->[0]->[1],14);
is($q_gapblocks->[1]->[1],21);
is($q_gapblocks->[1]->[0],1152);


is( $hit->next_hsp, undef, 'next_hsp should be undef');
is( $result->next_hit, undef, 'next_hit should be undef');
TODO: {
    local $TODO = "next_result should really return undef, not empty string";
    is( $pslparser->next_result, undef, 'next_result should be undef');
}

# bug 2850

my $searchio = Bio::SearchIO->new(
    -format => 'psl',
    -file   => test_input_file('headerless.psl'),
);

lives_ok { my $result = $searchio->next_result };
