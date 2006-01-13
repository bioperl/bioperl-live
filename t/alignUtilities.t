# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 9;

BEGIN {
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;

	plan tests => NUMTESTS;
}


use Bio::Align::Utilities qw(:all);
ok(1);
use Bio::SimpleAlign;
use Bio::PrimarySeq;
use Bio::LocatableSeq;



# hand crafting the simple input data
use Data::Dumper;

my $aa_align = new Bio::SimpleAlign;
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n1", -seq => "MLIDVG-MLVLR"));
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n1", -seq => "MLIDVRTPLALR"));
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n1", -seq => "MLI-VR-SLALR"));

my %dnaseqs = ();
$dnaseqs{'n1'} = Bio::PrimarySeq->new(-id => "n1", -seq => 'atgctgatagacgtaggcatgctagtactgagatga');
$dnaseqs{'n2'} = Bio::PrimarySeq->new(-id => "n2", -seq => 'atgctgatcgacgtacgcaccccgctagcactcagatga');
$dnaseqs{'n3'} = Bio::PrimarySeq->new(-id => "n3", -seq => 'atgttgattgtacgctcgcttgcacttagatga');

ok my $dna_aln = &aa_to_dna_aln($aa_align, \%dnaseqs);

#print Dumper $dna_aln;

ok $dna_aln->length, 36;
ok $dna_aln->no_residues, 99;
ok $dna_aln->no_sequences, 3;
ok $dna_aln->consensus_string(50), "atgctgatagacgtaggc????t??tacta?ta?ga";

$dna_aln->verbose(-1);
ok my $replicates = &bootstrap_replicates($dna_aln,3);

ok scalar @$replicates, 3;
my $repl_aln = pop @$replicates;
ok $repl_aln->no_sequences, 3;


#use IO::String;
#use Bio::AlignIO;
#my $string;
#my $out = IO::String->new($string);
#
#my $strout = Bio::AlignIO->new(-fh   => $out,'-format' => 'pfam');
#$strout->write_aln($repl_aln);
#ok $string, "";
