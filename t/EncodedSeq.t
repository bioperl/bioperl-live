
# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 40;

BEGIN {     
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
    
    plan tests => NUMTESTS;
	use_ok('Bio::Seq::EncodedSeq');
	use_ok('Bio::SimpleAlign');
	use_ok('Bio::AlignIO');
	use_ok('Bio::Root::IO');
}


my ($str, $aln, $seq, $loc);

ok $seq = new Bio::Seq::EncodedSeq(
			     -seq => '--atg---gta--',
			     -start => 1,
			     -end => 6,
			     -strand => 1
			     );
is $seq->alphabet, 'dna';
is $seq->start, 1;
is $seq->end, 6;
is $seq->strand, 1;
is $seq->no_gaps, 1;
is $seq->column_from_residue_number(4), 9;

# this should fail
eval {
    $seq->column_from_residue_number(8);
};
ok $@;

ok $loc = $seq->location_from_column(4);
isa_ok $loc, 'Bio::Location::Simple';
is $loc->to_FTstring, "2";

ok $loc = $seq->location_from_column(6);
isa_ok $loc,'Bio::Location::Simple';
is $loc->start, 3;
is $loc->location_type, 'IN-BETWEEN';
is $loc->to_FTstring, '3^4';

is $loc = $seq->location_from_column(2), undef;

is $seq->encoding, "GGCCCGGGCCCGG";
is $seq->encoding(-explicit => 1), "GGCDEGGGCDEGG";

ok $seq = new Bio::Seq::EncodedSeq(
			     -seq => 'atcgta',
			     -start => 10,
			     -end => 15,
			     -strand => -1,
			     );
is $seq->encoding('CCGGG'), 'CCGGGCCCC';
is $seq->seq, 'atcg---ta';
is $seq->column_from_residue_number(14), 2;
is $seq->encoding('3C2GCG'), 'CCCGGCGCC';
is $seq->seq, 'at-c--gta';
is $seq->no_gaps, 2;
is $seq->location_from_column(2)->to_FTstring, 14;
is $seq->location_from_column(5)->to_FTstring, "12^13";
is $seq->encoding("B", Bio::Location::Simple->new(-start => 10, -end => 11,
						  -location_type => 'IN-BETWEEN')), 'B';
is $seq->seq, 'at-c--gt-a';
is $seq->encoding, 'CBCCGGCGCC';
is $seq->cds(-nogaps => 1)->seq, 'tacgat';
is $seq->translate->seq, 'YD';
ok $seq = $seq->trunc(4,10); # kinda testing LocatableSeq's new trunc() here as well.
is $seq->seq, 'c--gt-a';
is $seq->encoding, 'CBCCGGC';
