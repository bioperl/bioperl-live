# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 22;

BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    
    plan tests => NUMTESTS;
}
use Bio::LocatableSeq;
ok(1);
use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::Root::IO;

my ($str, $aln, $seq, $loc);

ok $seq = new Bio::LocatableSeq(
			     -seq => '--atg---gta--',
			     -start => 1,
			     -end => 6,
			     -strand => 1,
			     -alphabet => 'dna'
			     );
ok $seq->alphabet, 'dna';
ok $seq->start, 1;
ok $seq->end, 6;
ok $seq->strand, 1;
ok $seq->no_gaps, 1;
ok $seq->column_from_residue_number(4), 9;

ok $loc = $seq->location_from_column(4);
ok $loc->isa('Bio::Location::Simple');
ok $loc->to_FTstring, 2;

ok $loc = $seq->location_from_column(6);
ok $loc->isa('Bio::Location::Simple');
ok $loc->start, 3;
ok $loc->location_type, 'IN-BETWEEN';
ok $loc->to_FTstring, '3^4';


ok $loc = $seq->location_from_column(2), undef;


$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $seq = $aln->get_seq_by_pos(1);
ok ref($seq), 'Bio::LocatableSeq';

ok $seq->get_nse, '1433_LYCES/9-246';
ok $seq->id, '1433_LYCES';

