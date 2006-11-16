# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
use strict;
use constant NUMTESTS => 13;
use vars qw($DEBUG);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;

BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;

	plan tests => NUMTESTS;
	use_ok('Bio::Align::Utilities', qw(:all));
	use_ok('Bio::SimpleAlign');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::LocatableSeq');
	use_ok('Bio::AlignIO');
}





# hand crafting the simple input data
use Data::Dumper;

my $aa_align = new Bio::SimpleAlign;
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n1", -seq => "MLIDVG-MLVLR"));
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n2", -seq => "MLIDVRTPLALR"));
$aa_align->add_seq(Bio::LocatableSeq->new(-id => "n3", -seq => "MLI-VR-SLALR"));

my %dnaseqs = ();
$dnaseqs{'n1'} = Bio::PrimarySeq->new(-id => "n1", -seq => 'atgctgatagacgtaggcatgctagtactgaga');
$dnaseqs{'n2'} = Bio::PrimarySeq->new(-id => "n2", -seq => 'atgctgatcgacgtacgcaccccgctagcactcaga');
$dnaseqs{'n3'} = Bio::PrimarySeq->new(-id => "n3", -seq => 'atgttgattgtacgctcgcttgcacttaga');
my $dna_aln;

ok( $dna_aln = &aa_to_dna_aln($aa_align, \%dnaseqs));
if( $DEBUG ) {
    Bio::AlignIO->new(-format=>'clustalw')->write_aln($dna_aln);
  }
#print Dumper $dna_aln;

is $dna_aln->length, 36;
is $dna_aln->no_residues, 99;
is $dna_aln->no_sequences, 3;
is $dna_aln->consensus_string(50), "atgctgat?gacgtacgc????cgctagcact?aga";

$dna_aln->verbose(-1);
my $replicates;
ok $replicates = &bootstrap_replicates($dna_aln,3);

is scalar @$replicates, 3;
my $repl_aln = pop @$replicates;
is $repl_aln->no_sequences, 3;


#use IO::String;
#use Bio::AlignIO;
#my $string;
#my $out = IO::String->new($string);
#
#my $strout = Bio::AlignIO->new(-fh   => $out,'-format' => 'pfam');
#$strout->write_aln($repl_aln);
#is $string, "";
