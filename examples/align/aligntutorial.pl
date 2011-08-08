#!/usr/bin/perl

# An example of how to use the different alignment tools in bioperl
# to align some sequences

# All these methods except Bio::Tools::pSW will work for DNA sequence
# (need to use a different matrix however)

use Bio::Factory::EMBOSS;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::pSW;
use Bio::PrimarySeq;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::TCoffee;
use Bio::Tools::Run::StandAloneBlast;

use strict;
# build the sequences since EMBOSS expects seqs to be in files
my $seq = new Bio::PrimarySeq(-seq => 
				  'MAVNPELAPFTLSRGIPSFDDQALSTIIQLQDCIQQAIQQLNYSTAEFLAELLYAECSILDKSSVYWSDAVYLYALSLFLNKSYHTAFQISKEFKEYHLGIAYIFGRCALQLSQGVNEAILTLLSIINVFSSNSSNTRINMVLNSNLVHIPDLATLNCLLGNLYMKLDHSKEGAFYHSEALAINPYLWESYEAICKMRATVDLKRVFFDIAGKKSNSHNNNAASSFPSTSLSHFEPRSQPSLYSKTNKNGNNNINNNVNTLFQSSNSPPSTSASSFSSIQHFSRSQQQQANTSIRTCQNKNTQTPKNPAINSKTSSALPNNISMNLVSPSSKQPTISSLAKVYNRNKLLTTPPSKLLNNDRNHQNNNNNNNNNNNNNNNNNNNNNNNNIINKTTFKTPRNLYSSTGRLTTSKKNPRSLIISNSILTSDYQITLPEIMYNFALILRSSSQYNSFKAIRLFESQIPSHIKDTMPWCLVQLGKLHFEIINYDMSLKYFNRLKDLQPARVKDMEIFSTLLWHLHDKVKSSNLANGLMDTMPNKPETWCCIGNLLSLQKDHDAAIKAFEKATQLDPNFAYAYTLQGHEHSSNDSSDSAKTCYRKALACDPQHYNAYYGLGTSAMKLGQYEEALLYFEKARSINPVNVVLICCCGGSLEKLGYKEKALQYYELACHLQPTSSLSKYKMGQLLYSMTRYNVALQTFEELVKLVPDDATAHYLLGQTYRIVGRKKDAIKELTVAMNLDPKGNQVIIDELQKCHMQE',
				  -id => 'seq1'
				  ); 
my $seq2 = new Bio::PrimarySeq( -seq =>
				'CLIFXRLLLIQMIHPQARRAFTFLQQQEPYRIQSMEQLSTLLWHLADLPALSHLSQSLISISRSSPQAWIAVGNCFSLQKDHDEAMRCFRRATQVDEGCAYAWTLCGYEAVEMEEYERAMAFYRTAIRTDARHYNAWYVLFFFFFFFFVPGDIDSXPKKGMEWGXFISKRIDRGMRSIILKEPSKSIQLIPFFYVALVWXVGVSSYPLETMTNIDFPKKKKALEKSNDVVQALHFYERASKYAPTSAMVQFKRIRALVALQRYDEAISALVPLTHSAPDEANVFFLLGKCLLKKERRQEATMAFTNARELEPK',
				    -id => 'seq2');

    my $out = new Bio::SeqIO(-format => 'fasta',
			     -file   => ">seq1.fa");
    $out->write_seq($seq);
    $out->close();
    $out = new Bio::SeqIO(-format => 'fasta',
			  -file   => ">seq2.fa");
    $out->write_seq($seq2);
    $out->close();


my $embossfactory = Bio::Factory::EMBOSS->new();

my @alignprogs = qw(water needle stretcher matcher);
my $alignout = new Bio::AlignIO(-format => 'msf');


foreach my $prog ( @alignprogs ) {
    my $alignfactory = $embossfactory->program('water');


    $alignfactory->run({ '-sequencea' => 'seq1.fa',
			    '-seqall'    => 'seq2.fa',
			    '-gapext'    => 2.0,
			    '-datafile'  => 'EBLOSUM62',
			    '-gapopen'   => 14.0,
			    '-outfile'   => "seq1_vs_seq2.$prog"});

    my $alnin = new Bio::AlignIO(-format => 'emboss',
				 -file   => "seq1_vs_seq2.$prog");
    my $aln = $alnin->next_aln(); 
    $alignout->write_aln($aln);
}

# this should produce the same alignment as 'water'
my $factory = new Bio::Tools::pSW(-matrix=> 'blosum62.bla',
				  -gap   => 14,
				  -ext   => 2);
my $aln = $factory->pairwise_alignment($seq,$seq2);
$alignout->write_aln($aln);

$factory = new Bio::Tools::Run::Alignment::Clustalw('ktuple' => 2, 
						    'matrix' => 'BLOSUM');
$aln = $factory->align([$seq,$seq2]);
$alignout->write_aln($aln);

$factory = new Bio::Tools::Run::Alignment::TCoffee('ktuple' => 2, 
						   'matrix' => 'BLOSUM');
$aln = $factory->align([$seq,$seq2]);
$alignout->write_aln($aln);

$factory = new Bio::Tools::Run::StandAloneBlast();
$aln = $factory->bl2seq($seq,$seq2);

# this actually returns a Bio::Tools::BPbl2seq object
# it can be transformed to a SimpleAlign object see
# the code in Bio::AlignIO::bl2seq
# A transformer object will be written at some point
