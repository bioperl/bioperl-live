# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 11;

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

my ($str, $aln, $seq);

ok $seq = new Bio::LocatableSeq(
			     -seq => '--atg---gta--',
			     -start => 1,
			     -end => 6
			     );

ok $seq->column_from_residue_number(4), 9;


exit;

$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $seq = $aln->get_seq_by_pos(1);
ok ref($seq), 'Bio::LocatableSeq';

ok $seq->get_nse, '1433_LYCES/9-246';
ok $seq->id, '1433_LYCES';
ok $seq->no_gaps, 1;


my $string = $seq->seq;
print $string, CORE::length($string), "\n";

use Data::Dumper;
print Dumper($seq);
