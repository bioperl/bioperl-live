
# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 31;

BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    
    plan tests => NUMTESTS;
}
use Bio::Seq::EncodedSeq;
ok(1);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::IO;

my ($str, $aln, $seq, $loc);

ok $seq = new Bio::Seq::EncodedSeq(
			     -seq => '--atg---gta--',
			     -start => 1,
			     -end => 6,
			     -strand => 1
			     );
ok $seq->alphabet, 'dna';
ok $seq->start, 1;
ok $seq->end, 6;
ok $seq->strand, 1;
ok $seq->no_gaps, 1;
ok $seq->column_from_residue_number(4), 9;
$@ = undef;
eval { $seq->column_from_residue_number(8) };
ok $@;

ok $loc = $seq->location_from_column(4);
ok $loc->isa('Bio::Location::Simple');
ok $loc->to_FTstring, 2;

ok $loc = $seq->location_from_column(6);
ok $loc->isa('Bio::Location::Simple');
ok $loc->start, 3;
ok $loc->location_type, 'IN-BETWEEN';
ok $loc->to_FTstring, '3^4';

ok $loc = $seq->location_from_column(2), undef;

ok $seq->encoding, "GGCCCGGGCCCGG";
ok $seq->encoding(-explicit => 1), "GGCDEGGGCDEGG";

ok $seq = new Bio::Seq::EncodedSeq(
			     -seq => 'atcgta',
			     -start => 10,
			     -end => 15,
			     -strand => -1,
			     );
ok $seq->encoding('CCGGG'), 'CCGGGCCCC';
ok $seq->seq, 'atcg---ta';
ok $seq->column_from_residue_number(14), 8;
ok $seq->encoding('3C2GCG'), 'CCCGGCGCC';
ok $seq->seq, 'at-c--gta';
ok $seq->no_gaps, 2;
ok $seq->location_from_column(2)->to_FTstring, 11;
ok $seq->location_from_column(5)->to_FTstring, "12^13";
ok $seq->encoding("B", Bio::Location::Simple->new(-start => 10, -end => 11,
						  -location_type => 'IN-BETWEEN')), 'B';
ok $seq->seq, 'a-t-c--gta';
ok $seq->encoding, 'CCCGGCGCBC';
ok $seq->cds(-nogaps => 1)->seq, 'tacgat';
ok $seq->translate->seq, 'YD';
ok $seq = $seq->trunc(4,10); # kinda testing LocatableSeq's new trunc() here as well.
ok $seq->seq, '-c--gta';
ok $seq->encoding, 'CCCGGCG';
