# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 45;

BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    
    plan tests => NUMTESTS;
}

use Bio::SimpleAlign;
ok(1);
use Bio::AlignIO;
use Bio::Root::IO;

my ($str, $aln, @seqs, $seq);

$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile("t","data","testaln.pfam"));
ok defined($str) && ref($str) && $str->isa('Bio::AlignIO');
$aln = $str->next_aln();
ok $aln->get_seq_by_pos(1)->get_nse, '1433_LYCES/9-246', " failed pfam input test";

my $aln2 = $aln->select(1,3);
ok $aln2;
ok $aln2->no_sequences, 3;

@seqs = $aln->each_seq();
ok scalar @seqs, 16;
ok $seqs[0]->get_nse, '1433_LYCES/9-246';
ok $seqs[0]->id, '1433_LYCES';
ok $seqs[0]->no_gaps, 3;
@seqs = $aln->each_alphabetically();
ok scalar @seqs, 16;

ok $aln->column_from_residue_number('1433_LYCES', 10), 2; 
ok $aln->displayname('1433_LYCES/9-246', 'my_seq'), 'my_seq';
ok $aln->displayname('1433_LYCES/9-246'), 'my_seq';
ok substr ($aln->consensus_string(50), 0, 60), 
    "RE??VY?AKLAEQAERYEEMV??MK?V????????ELS?EERNLLSVAYKNVIGARRASW";

ok (@seqs = $aln->each_seq_with_id('143T_HUMAN'));
ok scalar @seqs, 1;

ok $aln->is_flush, 1;
ok ($aln->id('x') and $aln->id eq 'x');

ok $aln->length, 242;
ok $aln->no_residues, 103;
ok $aln->no_sequences, 16;
ok $aln->percentage_identity(), 66.9052451661147 ;

ok $aln->set_displayname_count;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES_1';
ok $aln->set_displayname_flat;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES';
ok $aln->set_displayname_normal;
ok $aln->displayname('1433_LYCES/9-246'), '1433_LYCES/9-246';
ok $aln->uppercase;
ok $aln->map_chars('\.','-');
@seqs = $aln->each_seq_with_id('143T_HUMAN');
ok substr($seqs[0]->seq, 0, 60), 
    'KTELIQKAKLAEQAERYDDMATCMKAVTEQGA---ELSNEERNLLSVAYKNVVGGRRSAW';

ok $aln->remove_seq($seqs[0]);
ok $aln->no_sequences, 15;
ok $aln->add_seq($seqs[0]);
ok $aln->no_sequences, 16;
ok $seq = $aln->get_seq_by_pos(1);

ok (($aln->missing_char(), 'P') and  ($aln->missing_char('X'), 'X')) ;
ok (($aln->match_char(), '.') and  ($aln->match_char('-'), '-')) ;
ok (($aln->gap_char(), '-') and  ($aln->gap_char('.'), '.')) ;

# write test for:
# purge()

eval { require 'IO/String.pm' };
if( $@ ) {
    print STDERR "IO::String not installed.  Skipping tests.\n";
    for( 32..NUMTESTS ) {
	skip(1,"IO::String not installed. Skipping tests");
    }
    exit;
}


my $string;
my $out = IO::String->new($string);

my $s1 = new Bio::LocatableSeq (-id => 'AAA', 
			    -seq => 'aaaaattttt',
			    -start => 1,
			    -end => 10
			    );
my $s2 = new Bio::LocatableSeq (-id => 'BBB', 
			    -seq => '-aaaatttt-',
			    -start => 1,
			    -end => 8
			    );

$a = new Bio::SimpleAlign;
$a->add_seq($s1);
$a->add_seq($s2);


my $strout = Bio::AlignIO->new(-fh   => $out,'-format' => 'pfam');
$strout->write_aln($a);
ok $string, "AAA/1-10    aaaaattttt
BBB/1-8     -aaaatttt-
";

$out->setpos(0); $string ='';
my $b = $a->slice(2,9);
$strout->write_aln($b);
ok $string, "AAA/2-9    aaaatttt
BBB/1-8    aaaatttt
";

$out->setpos(0); $string ='';
$b = $a->slice(9,10);
$strout->write_aln($b);
ok $string, "AAA/9-10    tt
BBB/8-8     t-
";

$a->verbose(-1);
$out->setpos(0); $string ='';
$b = $a->slice(1,1);
$strout->write_aln($b);
ok $string, "AAA/1-1    a\n";

$out->setpos(0); $string ='';
$b = $a->slice(10,13);
$strout->write_aln($b);
ok $string, "AAA/10-10    t\n";

eval {
    $b = $a->slice(11,13);
};
ok 1 if $@;
