# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;

    plan tests => 32;
}

use Bio::SimpleAlign;
ok(1);
use Bio::AlignIO;
use Bio::Root::IO;

my ($str, $aln, @seqs);

$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->{order}->{'0'}, '1433_LYCES/9-246', " failed pfam input test";

#use Data::Dumper;
#print Dumper($aln);

@seqs = $aln->each_seq();
ok scalar @seqs, 16;
ok $seqs[0]->get_nse, '1433_LYCES/9-246';
ok $seqs[0]->id, '1433_LYCES';
@seqs = $aln->each_alphabetically();
ok scalar @seqs, 16;

ok $aln->column_from_residue_number('1433_LYCES', 10), 2; 
ok $aln->displayname('1433_LYCES/9-246', 'my_seq'), 'my_seq';
ok $aln->displayname('1433_LYCES/9-246'), 'my_seq';
ok $aln->consensus_string(50), 
    "RE??VY?AKLAEQAERYEEMV??MK?V????????ELS?EERNLLSVAYKNVIGARRASWRIISSIEQK";

ok (@seqs = $aln->each_seq_with_id('143T_HUMAN'));
ok scalar @seqs, 1;

ok $aln->is_flush, 1;
ok ($aln->id('x') and $aln->id eq 'x');

ok $aln->length, 69;
ok $aln->no_residues, 48;
ok $aln->no_sequences, 16;
ok $aln->percentage_identity(), 71.2440864339599 ;

ok $aln->set_displayname_count;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES_1';
ok $aln->set_displayname_flat;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES';
ok $aln->set_displayname_normal;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES/9-246';
ok $aln->uppercase;
ok $aln->map_chars('\.','-');
@seqs = $aln->each_seq_with_id('143T_HUMAN');
ok $seqs[0]->seq, 'KTELIQKAKLAEQAERYDDMATCMKAVTEQGA---ELSNEERNLLSVAYKNVVGGRRSAWRVISSIEQK';

ok $aln->remove_seq($seqs[0]);
ok $aln->no_sequences, 15;
ok $aln->add_seq($seqs[0]);
ok $aln->no_sequences, 16;

# write test for:
# purge()
