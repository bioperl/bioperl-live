#!/usr/bin/env perl
use rlib '.';
use strict; use warnings;
use Test::More;
use Config;
use Helper;

my %notes = (
    # 'random' => 'build tree of random subset',
    'dist-all'    => 'half-matrix list of distances between leaves',
    'label-nodes' => 'prepends ID to leaf/node labels',
    'length'      => 'total branch length',
    'length-all'  => 'nodes and branch lengths',
    'otus-num'    => 'total number of OTUs',
    'otus-all'    => 'leaf nodes with branch lengths',
    # 'rmbl' => 'remove branch lengths from tree',
);

test_no_arg_opts('biotree', 'test-biotree.dnd', \%notes);

my $opts = [
    ['otus-desc', '15',
     "Prints OTU's descended from node 15"],
    ['depth', 'N40,B31,SV1', 'depth to root'],
    ['dist', 'SV1,N40',
     'distance between a pair of nodes or leaves'],
    ['lca', 'SV1,B31,N40', 'Least Common Ancestor'],
    ['ltt', '10',
     'divide tree into 10 segments and count branches'],
    ['output', 'tabtree', 'Output file format tabtree'],
    ['subset', 'SV1,B31,N40', 'subset tree'],
    ['walk', '156a', 'walk tree from 156a'],
    ];

test_one_arg_opts('biotree', 'test-biotree.dnd', $opts);

run_bio_program_nocheck('biotree', 'test-biotree.dnd', "--random 1",
			{note=>"Random of sample size 1"});

note( "Skipping biotree -r JD1 - needs investigation" );
# run_bio_program_nocheck('biotree', 'test-biotree.dnd', "-r JD1",
#			{note=>"resample"});

done_testing();
