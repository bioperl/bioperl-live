# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 66;

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

my $aln1 = $aln->remove_columns(['mismatch']);
ok ($aln1->match_line, '::*::::*:**:*:*:***:**.***::*.*::**::**:***..**:*:*.::::*:.:*.*.**:***.**:*.:.**::**.*:***********:::*:.:*:**.*::*:.*.:*:**:****************::');

my $aln2 = $aln->select(1,3);
ok $aln2;
ok $aln2->no_sequences, 3;

# test select non continous
$aln2 = $aln->select_noncont(8,2,7);
ok($aln2->no_sequences, 3);
ok($aln2->get_seq_by_pos(2)->id, $aln->get_seq_by_pos(7)->id);

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
    "RE??VY?AKLAEQAERYEEMV??MK?VAE??????ELSVEERNLLSVAYKNVIGARRASW";
ok substr ($aln->consensus_string(100), 0, 60),
    "?????????L????E????M???M????????????L??E?RNL?SV?YKN??G??R??W";
ok substr ($aln->consensus_string(0), 0, 60),
    "REDLVYLAKLAEQAERYEEMVEFMKKVAELGAPAEELSVEERNLLSVAYKNVIGARRASW";

ok (@seqs = $aln->each_seq_with_id('143T_HUMAN'));
ok scalar @seqs, 1;

ok $aln->is_flush, 1;
ok ($aln->id('x') and $aln->id eq 'x');

ok $aln->length, 242;
ok $aln->no_residues, 3769;
ok $aln->no_sequences, 16;
ok (sprintf("%.2f",$aln->overall_percentage_identity()), 33.06);
ok (sprintf("%.2f",$aln->average_percentage_identity()), 66.91);

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

ok($aln->match_line, '       ::*::::*  : *   *:           *: *:***:**.***::*. *::**::**:***      .  .      **  :* :*   .  :: ::   *:  .     :* .*. **:***.** :*.            :  .*  *   :   : **.*:***********:::* : .: *  :** .*::*: .*. : *: **:****************::     ');
ok $aln->remove_seq($seqs[0]);
ok $aln->no_sequences, 15;
ok $aln->add_seq($seqs[0]);
ok $aln->no_sequences, 16;
ok $seq = $aln->get_seq_by_pos(1);
ok( $seq->id, '1433_LYCES');
ok (($aln->missing_char(), 'P') and  ($aln->missing_char('X'), 'X')) ;
ok (($aln->match_char(), '.') and  ($aln->match_char('-'), '-')) ;
ok (($aln->gap_char(), '-') and  ($aln->gap_char('.'), '.')) ;

ok $aln->purge(0.7), 12;
ok $aln->no_sequences, 4;

eval { require 'IO/String.pm' };
if( $@ ) {
	print STDERR "IO::String not installed.  Skipping tests.\n";
	for( $Test::ntest..NUMTESTS ) {
		skip("IO::String not installed. Skipping tests",1);
	}
	exit;
}

my $string;
my $out = IO::String->new($string);

my $s1 = new Bio::LocatableSeq (-id => 'AAA',
			    -seq => 'aawtat-tn-',
			    -start => 1,
			    -end => 8,
  			    -alphabet => 'dna'
			    );
my $s2 = new Bio::LocatableSeq (-id => 'BBB',
			    -seq => '-aaaat-tt-',
			    -start => 1,
			    -end => 7,
  			    -alphabet => 'dna'
			    );
$a = new Bio::SimpleAlign;
$a->add_seq($s1);           
$a->add_seq($s2);

ok $a->consensus_iupac, "aAWWAT-TN-";
$s1->seq('aaaaattttt');
$s1->alphabet('dna');
$s1->end(10);
$s2->seq('-aaaatttt-');
$s2->end(8);
$a = new Bio::SimpleAlign;
$a->add_seq($s1);
$a->add_seq($s2);

my $strout = Bio::AlignIO->new(-fh   => $out,'-format' => 'pfam');
$strout->write_aln($a);
ok $string, "AAA/1-10    aaaaattttt
BBB/1-8     -aaaatttt-
";

$out->setpos(0); 
$string ='';
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
$b = $a->slice(1,2);
$strout->write_aln($b);
ok $string, "AAA/1-2    aa
BBB/1-1    -a
";

eval {
	$b = $a->slice(11,13);
};

ok ($@ =~ /EX/ );

# remove_columns by position
$out->setpos(0); 
$string ='';
$str = Bio::AlignIO->new(-file=> Bio::Root::IO->catfile(
											"t","data","mini-align.aln"));
$aln1 = $str->next_aln;
$aln2 = $aln1->remove_columns([0,0]);
$strout->write_aln($aln2);
ok $string, "P84139/1-33              NEGEHQIKLDELFEKLLRARLIFKNKDVLRRC
P814153/1-33             NEGMHQIKLDVLFEKLLRARLIFKNKDVLRRC
BAB68554/1-14            ------------------AMLIFKDKQLLQQC
gb|443893|124775/1-32    MRFRFQIKVPPAVEGARPALLIFKSRPELGGC
";

# and when arguments are entered in "wrong order"?
$out->setpos(0); 
$string ='';
my $aln3 = $aln1->remove_columns([1,1],[30,30],[5,6]);
$strout->write_aln($aln3);
ok $string, "P84139/1-33              MEGEIKLDELFEKLLRARLIFKNKDVLRC
P814153/1-33             MEGMIKLDVLFEKLLRARLIFKNKDVLRC
BAB68554/1-14            ----------------AMLIFKDKQLLQC
gb|443893|124775/1-32    -RFRIKVPPAVEGARPALLIFKSRPELGC
";

my %cigars = $aln1->cigar_line;
ok $cigars{'gb|443893|124775/1-32'},'19,19:21,24:29,29:32,32';
ok $cigars{'P814153/1-33'},'20,20:22,25:30,30:33,33';
ok $cigars{'BAB68554/1-14'},'1,1:3,6:11,11:14,14';
ok $cigars{'P84139/1-33'},'20,20:22,25:30,30:33,33';


# sort_alphabetically
my $s3 = new Bio::LocatableSeq (-id => 'ABB',
										  -seq => '-attat-tt-',
										  -start => 1,
										  -end => 7,
										  -alphabet => 'dna'
										 );
$a->add_seq($s3);

ok $a->get_seq_by_pos(2)->id,"BBB";
ok $a->sort_alphabetically;
ok $a->get_seq_by_pos(2)->id,"ABB";

$b = $a->remove_gaps();
ok $b->consensus_string, "aaaattt";

$s1->seq('aaaaattt--');

$b = $a->remove_gaps(undef, 'all_gaps_only');
ok $b->consensus_string, "aaaaatttt";

__END__

  print $aln->score;
  print $aln->percentage_identity;
