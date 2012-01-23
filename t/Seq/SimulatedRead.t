# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use warnings;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    test_begin(-tests => 194);

    use_ok('Bio::Seq');
    use_ok('Bio::Seq::Quality');
    use_ok('Bio::PrimarySeq');
    use_ok('Bio::LocatableSeq');
    use_ok('Bio::Seq::SimulatedRead');
}

my $VERBOSE = test_debug();

my ($ref, $ref2, $ref3, $ref4, $ref5, $read, $errors);

$ref = Bio::Seq::Quality->new(-id    => 'human_id',
                               -seq   => 'TAAAAAAACCCC',
                               -qual  => '1 2 3 4 5 6 7 8 9 10 11 12',
                               -trace => '0 5 10 15 20 25 30 35 40 45 50 55',
                               -desc  => 'The human genome' );

$ref2 = Bio::Seq->new(-id   => 'other_genome',
                       -seq  => 'ACGTACGT',
                       -desc => '"Secret" sequence');

$ref3 = Bio::PrimarySeq->new(-seq => 'ACACTGATCTAGCGTCGTGCTAGCTGACGTAGCTGAT' );

$ref4 = Bio::LocatableSeq->new(-id  => 'a_thaliana',
                                -seq => 'CGTATTCTGAGGAGAGCTCT' );


# Basic object

ok $read = Bio::Seq::SimulatedRead->new();
isa_ok $read, 'Bio::Seq::SimulatedRead';
isa_ok $read, 'Bio::LocatableSeq';
isa_ok $read, 'Bio::Seq::Quality';

$errors->{'1'}->{'+'} = 'T';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->reference, $ref;
ok $read->errors;
is $read->errors->{'1'}->{'+'}->[0], 'T';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -track => 0 );
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is $read->track, 0;
is $read->desc, undef;
is $read->revcom->seq, 'GGGGTTTTTTTA';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -track => 1 );
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is join(' ',@{$read->qual}), '';
is $read->track, 1;
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -track => 1, -coord_style => 'bioperl' );
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -track => 1, -coord_style => 'genbank' );
is $read->desc, 'reference=human_id position=1..12 description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -qual_levels => [30, 10]);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30 30 30 30 30 30';
is $read->track, 1;
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';
is $read->revcom->seq, 'GGGGTTTTTTTA';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref2 );
is $read->start, 1;
is $read->end, 8;
is $read->seq, 'ACGTACGT';
is join(' ',@{$read->qual}), '';
is $read->desc, 'reference=other_genome start=1 end=8 strand=+1 description="\"Secret\" sequence"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref3 );
is $read->start, 1;
is $read->end, 37;
is $read->seq, 'ACACTGATCTAGCGTCGTGCTAGCTGACGTAGCTGAT';
is join(' ',@{$read->qual}), '';
is $read->desc, 'start=1 end=37 strand=+1';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => -1, -qual_levels => [30, 10]);
is $read->seq, 'GGGGTTTTTTTA';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=-1 description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => -1, -qual_levels => [30, 10], -coord_style => 'genbank' );
is $read->desc, 'reference=human_id position=complement(1..12) description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -start => 2, -end => 8, -qual_levels => [30, 10]);
is $read->seq, 'AAAAAAA';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=2 end=8 strand=+1 description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => -1, -start => 2, -end => 8, -qual_levels => [30, 10]);
is $read->seq, 'TTTTTTT';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=2 end=8 strand=-1 description="The human genome"';

$errors = {};
$errors->{'6'}->{'+'} = 'GG';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => -1, -start => 2, -end => 8, -errors => $errors, -qual_levels => [30, 10]);
is $read->start, 2;
is $read->end, 8;
is $read->seq, 'TTTTTTGGT';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 10 10 30';
is $read->desc, 'reference=human_id start=2 end=8 strand=-1 errors=6+G,6+G description="The human genome"';

$errors = {};
$errors->{'6'}->{'+'} = 'GG';
$errors->{'1'}->{'%'} = 'T';
$errors->{'3'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => 1, -start => 2, -end => 8, -errors => $errors, -qual_levels => [30, 10]);
is $read->start, 2;
is $read->end, 8;
is $read->seq, 'TAAAAGGA';
is join(' ', @{$read->qual}), '10 30 30 30 30 10 10 30';
is $read->desc, 'reference=human_id start=2 end=8 strand=+1 errors=1%T,3-,6+G,6+G description="The human genome"';

$errors = {};
$errors->{'6'}->{'+'} = 'GG';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors, -qual_levels => [30, 10]);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAGGAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 10 10 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 errors=6+G,6+G description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors, -mid => 'ACGT', -errors => $errors, -qual_levels => [30, 10]);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'ACGTTAGGAAAAAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 10 10 30 30 30 30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 mid=ACGT errors=6+G,6+G description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -mid => 'TTTAAA', -qual_levels => [30, 10]);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TTTAAATAAAAAAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 mid=TTTAAA description="The human genome"';

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -mid => '', -qual_levels => []);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is join(' ', @{$read->qual}), '';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';


# Redundant errors

$errors = {};
$errors->{'6'}->{'+'} = ['G', 'G'];
$errors->{'1'}->{'%'} = ['A', 'G', 'T'];
$errors->{'3'}->{'-'} = [undef, undef];
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -strand => 1, -start => 2, -end => 8, -errors => $errors, -qual_levels => [30, 10]), 'redundant errors';
is $read->start, 2;
is $read->end, 8;
is $read->seq, 'TAAAAGGA';
is join(' ', @{$read->qual}), '10 30 30 30 30 10 10 30';
is $read->desc, 'reference=human_id start=2 end=8 strand=+1 errors=1%A,1%G,1%T,3-,3-,6+G,6+G description="The human genome"';


# Specifying errors() after new()

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -mid => '', -qual_levels => []);
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is join(' ', @{$read->qual}), '';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';

$errors = {};
ok $read->errors($errors), 'errors()';
is $read->start, 1;
is $read->end, 12;
is $read->seq, 'TAAAAAAACCCC';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 description="The human genome"';

$errors = {};
$errors->{'6'}->{'+'} = 'GG';
ok $read->errors($errors);
is $read->seq, 'TAAAAAGGAACCCC';
is $read->start, 1;
is $read->end, 12;
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 errors=6+G,6+G description="The human genome"';


# More tracking tests

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -mid => 'ACGT', -qual_levels => [], -coord_style => 'genbank' );
is $read->desc, 'reference=human_id position=1..12 mid=ACGT description="The human genome"';

ok $read->mid('AAAA');
is $read->desc, 'reference=human_id position=1..12 mid=AAAA description="The human genome"';

$errors = {};
$errors->{'6'}->{'+'} = 'GG';
ok $read->errors($errors);
is $read->desc, 'reference=human_id position=1..12 mid=AAAA errors=6+G,6+G description="The human genome"';

ok not($read->track(0)), 'track()';
is $read->track, 0;
is $read->desc, undef;
ok $read->track(1);
is $read->track, 1;
is $read->desc, 'reference=human_id position=1..12 mid=AAAA errors=6+G,6+G description="The human genome"';


# qual_levels() method

ok $read = Bio::Seq::SimulatedRead->new(-verbose => $VERBOSE, );
ok $read->qual_levels([30, 10]), 'qual_levels()';
is join(' ', @{$read->qual_levels}), '30 10';

# reference() method

ok $read->reference($ref), 'reference()';
is $read->reference(), $ref;

# mid() method

ok $read = Bio::Seq::SimulatedRead->new(-verbose => $VERBOSE, ), 'mid()';
ok $read->qual_levels([30, 10]);
ok $read->reference($ref);
ok $read->mid('ACGT');
ok $read->mid, 'ACGT';

is $read->seq, 'ACGTTAAAAAAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 mid=ACGT description="The human genome"';

ok $read->mid('TTTAAA');
ok $read->mid, 'TTTAAA';
is $read->seq, 'TTTAAATAAAAAAACCCC';
is join(' ', @{$read->qual}), '30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30';
is $read->desc, 'reference=human_id start=1 end=12 strand=+1 mid=TTTAAA description="The human genome"';


# Edge case... mutation of the last bases of a simulated read with MID

$errors = {};
$errors->{'18'}->{'%'} = 'T';
$read->errors($errors);
is $read->seq, 'TTTAAATAAAAAAACCCT';


# Try different BioPerl object types

ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref  ), 'Bio::Seq::Quality';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref2 ), 'Bio::Seq';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref3 ), 'Bio::PrimarySeq';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref4 ), 'Bio::LocatableSeq';


# More detailed tests of the error specifications

$errors = {};
$errors->{'0'}->{'-'} = undef;
warning_like {$read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors )}
    qr/Positions of substitutions and deletions have to be strictly positive but got 0/;
is $read->seq, 'TAAAAAAACCCC';

$errors = {};
$errors->{'1'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'AAAAAAACCCC';

$errors = {};
$errors->{'12'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TAAAAAAACCC';

$errors = {};
$errors->{'13'}->{'-'} = undef;
warning_like {$read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors )}
    qr/Position 13 is beyond end of read \(12 residues\)/; # there should be a warning too
is $read->seq, 'TAAAAAAACCCC';

$errors = {};
$errors->{'0'}->{'%'} = 'G';
warning_like {$read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors )}
    qr/Positions of substitutions and deletions have to be strictly positive/; # there should be a warning too
is $read->seq, 'TAAAAAAACCCC';

$errors = {};
$errors->{'1'}->{'%'} = 'G';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'GAAAAAAACCCC';

$errors = {};
$errors->{'12'}->{'%'} = 'G';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TAAAAAAACCCG';

$errors = {};
$errors->{'13'}->{'%'} = 'G';
warning_like { $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors ) }
    qr/Position 13 is beyond end of read \(12 residues\)/; # there should be a warning too
is $read->seq, 'TAAAAAAACCCC';

$errors = {};
$errors->{'0'}->{'+'} = 'A';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'ATAAAAAAACCCC';

$errors = {};
$errors->{'1'}->{'+'} = 'A';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TAAAAAAAACCCC';

$errors = {};
$errors->{'12'}->{'+'} = 'A';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TAAAAAAACCCCA';

$errors = {};
$errors->{'13'}->{'+'} = 'A';
warning_like {$read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors ) }
    qr/Position 13 is beyond end of read \(12 residues\)/; # there should be a warning too
is $read->seq, 'TAAAAAAACCCC';

$errors = {};
$errors->{'1'}->{'%'} = 'G';
$errors->{'2'}->{'%'} = 'G';
$errors->{'3'}->{'%'} = 'G';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'GGGAAAAACCCC';

$errors = {};
$errors->{'1'}->{'+'} = 'G';
$errors->{'2'}->{'+'} = 'G';
$errors->{'3'}->{'+'} = 'G';
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TGAGAGAAAAACCCC';

$errors = {};
$errors->{'1'}->{'-'} = undef;
$errors->{'2'}->{'-'} = undef;
$errors->{'3'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'AAAAACCCC';

$errors = {};
$errors->{'1'}->{'+'} = 'GGG';
$errors->{'2'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TGGGAAAAAACCCC';

$errors = {};
$errors->{'2'}->{'+'} = 'CC';
$errors->{'2'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TCCAAAAAACCCC';

$errors = {};
$errors->{'2'}->{'%'} = 'C';
$errors->{'2'}->{'-'} = undef;
ok $read = Bio::Seq::SimulatedRead->new(-reference => $ref, -errors => $errors );
is $read->seq, 'TAAAAAACCCC';
