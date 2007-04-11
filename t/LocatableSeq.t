# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 84;

BEGIN {
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => NUMTESTS;
	use_ok('Bio::LocatableSeq');
	use_ok('Bio::SimpleAlign');
	use_ok('Bio::AlignIO');
	use_ok('Bio::Root::IO');
}

my ($str, $aln, $seq, $loc);

ok $seq = new Bio::LocatableSeq(
			     -seq => '--atg---gta--',
			     -strand => 1,
			     -alphabet => 'dna'
			     );
is $seq->alphabet, 'dna';
is $seq->start, 1;
is $seq->end, 6;
is $seq->strand, 1;
is $seq->no_gaps, 1;
is $seq->column_from_residue_number(4), 9;

ok $loc = $seq->location_from_column(4);
isa_ok $loc,'Bio::Location::Simple';
is $loc->to_FTstring, 2;

ok $loc = $seq->location_from_column(6);
isa_ok $loc,'Bio::Location::Simple';
is $loc->start, 3;
is $loc->location_type, 'IN-BETWEEN';
is $loc->to_FTstring, '3^4';


is $loc = $seq->location_from_column(2), undef;


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str);
isa_ok $str,'Bio::AlignIO';
$aln = $str->next_aln();
ok $seq = $aln->get_seq_by_pos(1);
is ref($seq), 'Bio::LocatableSeq';

is $seq->get_nse, '1433_LYCES/9-246';
is $seq->id, '1433_LYCES';

# test revcom and trunc

$seq = new Bio::LocatableSeq(
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
#use Data::Dumper;
#print Dumper $seq;
#print Dumper $seq2;
#exit;
$seq2 = $seq->revcom();
is $seq2->seq, '--tac---cat--';
is $seq2->start, $seq->start;
is $seq2->end, $seq->end;
is $seq2->strand, $seq->strand * -1;

# test column-mapping for -1 strand sequence
$seq = new Bio::LocatableSeq(
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


ok $seq = new Bio::LocatableSeq(
			     -seq => '--atg---gta--',
			     -strand => -1,
			     -alphabet => 'dna'
			     );
is $seq->alphabet, 'dna';
is $seq->start, 1;
is $seq->end, 6;
is $seq->strand, -1;
is $seq->no_gaps, 1;
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

ok $seq = new Bio::LocatableSeq();
is $seq->seq, undef;
is $seq->start, undef;
is $seq->end, undef;
