# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 119);

    use_ok('Bio::LocatableSeq');
    use_ok('Bio::AlignIO');
}

my ($str, $aln, $seq, $loc);

# basic tests

ok $seq = Bio::LocatableSeq->new(
                 -seq => '--atg---gta--',
                 -strand => 1,
                 -alphabet => 'dna'
                 );
is $seq->alphabet, 'dna';
is $seq->start, 1;
is $seq->end, 6;
is $seq->strand, 1;
is $seq->num_gaps, 1;
is $seq->column_from_residue_number(4), 9;
is $seq->column_from_residue_number(3), 5;

ok $loc = $seq->location_from_column(4);
isa_ok $loc,'Bio::Location::Simple';
is $loc->to_FTstring, 2;

ok $loc = $seq->location_from_column(6);
isa_ok $loc,'Bio::Location::Simple';
is $loc->start, 3;
is $loc->location_type, 'IN-BETWEEN';
is $loc->to_FTstring, '3^4';

is $loc = $seq->location_from_column(2), undef;
TODO: {
  local $TODO = "Need to fix columns before start of seq w/ start > 1";
  $seq->start(90);
  is $loc = $seq->location_from_column(2), undef;
}

$str = Bio::AlignIO->new(-file=> test_input_file('testaln.pfam'));
ok defined($str);
isa_ok $str,'Bio::AlignIO';
$aln = $str->next_aln();
ok $seq = $aln->get_seq_by_pos(1);
is ref($seq), 'Bio::LocatableSeq';

is $seq->get_nse, '1433_LYCES/9-246';
is $seq->id, '1433_LYCES';

# test invalid sequence

throws_ok{ $seq = Bio::LocatableSeq->new( -seq => '//!\\' ) } qr/.+/;

# test revcom and trunc

$seq = Bio::LocatableSeq->new(
                 -seq => '--atg---gta--',
                 -strand => 1,
                 -alphabet => 'dna'
                 );

my $seq2 = $seq->trunc(1,9);
is $seq2->seq, '--atg---g';
is $seq2->start, 1;
is $seq2->end, 4;
is $seq2->strand, $seq->strand;

$seq2 = $seq->trunc(3,8);
is $seq2->seq, 'atg---';
is $seq2->start, 1;
is $seq2->end, 3;

is $seq->strand(-1), -1;
is $seq->start, 1;
is $seq->end, 6;
$seq2 = $seq->trunc(3,8);
is $seq2->seq, 'atg---';
is $seq2->start, 4;
is $seq2->end, 6;
$seq2 = $seq->revcom();
is $seq2->seq, '--tac---cat--';
is $seq2->start, $seq->start;
is $seq2->end, $seq->end;
is $seq2->strand, $seq->strand * -1;
is $seq2->column_from_residue_number(4), 9;
is $seq2->column_from_residue_number(3), 5;

# test column-mapping for -1 strand sequence
$seq = Bio::LocatableSeq->new(
                 -seq => '--atg---gtaa-',
                 -strand => -1,
                 -alphabet => 'dna'
                 );
is $seq->column_from_residue_number(5),5;
is $seq->column_from_residue_number(4),9;
ok $loc = $seq->location_from_column(4);
isa_ok $loc,'Bio::Location::Simple';
is $loc->to_FTstring, 6;
ok $loc = $seq->location_from_column(6);
isa_ok $loc,'Bio::Location::Simple';
is $loc->start, 4;
is $loc->location_type, 'IN-BETWEEN';
is $loc->to_FTstring, '4^5';


# more tests for trunc() with strand -1


ok $seq = Bio::LocatableSeq->new(
                 -seq => '--atg---gta--',
                 -strand => -1,
                 -alphabet => 'dna'
                 );
is $seq->alphabet, 'dna';
is $seq->start, 1;
is $seq->end, 6;
is $seq->strand, -1;
is $seq->num_gaps, 1;
is $seq->column_from_residue_number(4), 5;


ok $seq2 = $seq->trunc(1,9);
is $seq2->seq, '--atg---g';
is $seq2->start, 3;
is $seq2->end, 6;
is $seq2->strand, $seq->strand;

is $seq->location_from_column(3)->start, 6;
is $seq->location_from_column(11)->start, 1;
is $seq->location_from_column(9)->start, 3;



ok $seq2 = $seq->trunc(7,12);
is $seq2->seq, '--gta-';
is $seq2->start, 1;
is $seq2->end, 3;


ok $seq2 = $seq->trunc(2,6);
is $seq2->seq, '-atg-';
is $seq2->start, 4;
is $seq2->end, 6;

ok $seq2 = $seq->trunc(4,7);
is $seq2->seq, 'tg--';
is $seq2->start, 4;
is $seq2->end, 5;

ok $seq = Bio::LocatableSeq->new();
is $seq->seq, undef;
is $seq->start, undef;
is $seq->end, undef;
my $nse;
eval{$nse = $seq->get_nse};
ok($@);
is ($nse, undef);
$seq->force_nse(1);
eval{$nse = $seq->get_nse};
ok(!$@);
is ($nse, '/0-0');

# test mapping

# mapping only supported for 1 => 1, 3 => 1, or 1 => 3 mapping relationships

eval{$seq = Bio::LocatableSeq->new(
                 -mapping => [40 => 2],
                 );};

ok($@);
like($@, qr/Mapping values other than 1 or 3 are not currently supported/);

eval{$seq = Bio::LocatableSeq->new(
                 -mapping => [3 => 3],
                 );};

ok($@);

# sequence is translated to protein, retains original DNA coordinates
# mapping is 1 residue for every 3 coordinate positions
$seq = Bio::LocatableSeq->new(
                 -seq => 'KKKAIDLVGVDKARENRQAIYLGASAIAEF',
                 -strand => -1,
                 -mapping => [1 => 3],
                 -start => 1,
                 -end => 90,
                 -alphabet => 'dna'
                 );

is $seq->seq, 'KKKAIDLVGVDKARENRQAIYLGASAIAEF';
is $seq->start, 1;
is $seq->end, 90;

# sequence is reverse-translated to DNA, retains original protein coordinates
# mapping is 3 residues for every 1 coordinate positions
$seq = Bio::LocatableSeq->new(
                 -seq => 'aaraaraargcnathgayytngtnggngtngayaargcnmgngaraaymgncargcnathtayytnggngcnwsngcnathgcngartty',
                 -strand => -1,
                 -mapping => [3 => 1],
                 -start => 1,
                 -end => 30,
                 -alphabet => 'protein'
                 );

is $seq->seq, 'aaraaraargcnathgayytngtnggngtngayaargcnmgngaraaymgncargcnathtayytnggngcnwsngcnathgcngartty';
is $seq->start, 1;
is $seq->end, 30;

# frameshifts (FASTA-like)
# support for this is preliminary
# this is a real example from a TFASTY report

$seq = Bio::LocatableSeq->new(
                 -seq => 'MGSSSTDRELLSAADVGRTVSRIAHQIIEKTALDDPAERTRVVLLGIPTRGVILATRLAAKIKEFAGEDVPHGALDITLYRDDLNFKPPRPLEATSIPAF\GGVDDAIVILVDDVLYSGRSVRSALDALRDIGRPRIVQLAVLVDRGHRELPI--/DYVGKNVPTSRSESVHVLLSEHDDRDGVVISK',
                 -strand => 1,
                 -mapping => [1 => 3],
                 -start => 1,
                 -end => 552,
                 -frameshifts => { # position, frameshift
                    298 => -1,
                    455 => 1
                    },
                 -alphabet => 'dna'
                 );

is $seq->seq, 'MGSSSTDRELLSAADVGRTVSRIAHQIIEKTALDDPAERTRVVLLGIPTRGVILATRLAAKIKEFAGEDVPHGALDITLYRDDLNFKPPRPLEATSIPAF\GGVDDAIVILVDDVLYSGRSVRSALDALRDIGRPRIVQLAVLVDRGHRELPI--/DYVGKNVPTSRSESVHVLLSEHDDRDGVVISK';
is $seq->start, 1;
is $seq->end, 552;
$seq->verbose(2);
eval { $seq->end(554);};
ok $@;
like $@, qr/Overriding value \[554\] with value 552/;

lives_ok { $seq = Bio::LocatableSeq->new(
                 -seq => 'LSYC*',
                 -strand => 0,
                 -start => 1,
                 -end => 5,
                 -verbose => 2
                 );} '* is counted in length';

throws_ok { $seq = Bio::LocatableSeq->new(
                 -seq => 'LSYC*',
                 -strand => 0,
                 -start => 1,
                 -end => 6,
                 -verbose => 2
                 );} qr/Overriding value \[6\] with value 5/, '* is counted in length, but end is wrong';

# setting symbols (class variables) - demonstrate scoping issues when using
# globals with and w/o localization.  To be fixed in a future BioPerl version

# see bug 2715
my $temp;

{
    $temp = $Bio::LocatableSeq::GAP_SYMBOLS;
    $Bio::LocatableSeq::GAP_SYMBOLS = '-\?';
    $seq = Bio::LocatableSeq->new(
                     -seq => '??atg-?-gta-?',
                     -strand => 1,
                     -start => 10,
                     -end => 15,
                     -alphabet => 'dna',
                     );
    is $Bio::LocatableSeq::GAP_SYMBOLS, '-\?';    
    is $seq->start, 10;
    is $seq->end, 15;
}

is $Bio::LocatableSeq::GAP_SYMBOLS, '-\?';
is $seq->end(15), 15;
$Bio::LocatableSeq::GAP_SYMBOLS = $temp;
is $Bio::LocatableSeq::GAP_SYMBOLS, '\-\.=~';

{
    local $Bio::LocatableSeq::GAP_SYMBOLS = '-\?';
    $seq = Bio::LocatableSeq->new(
                     -seq => '??atg-?-gta-?',
                     -strand => 1,
                     -start => 10,
                     -end => 15,
                     -alphabet => 'dna',
                     );
    is $Bio::LocatableSeq::GAP_SYMBOLS, '-\?';    
    is $seq->start, 10;
    is $seq->end, 15;
}

is $seq->end, 15;

# note, recalling the end() method uses old $GAP_SYMBOLS, which
# no longer are set (this argues for locally set symbols)
TODO: {
    local $TODO = 'Bio::LocatableSeq global variables have scoping issues';
    is $Bio::LocatableSeq::GAP_SYMBOLS, '-\?';
    # this should be 15 
    isnt $seq->end(19), 19;
}

