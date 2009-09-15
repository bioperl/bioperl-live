# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13);
	
	use_ok('Bio::Align::Utilities', qw(:all));
	use_ok('Bio::SimpleAlign');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::LocatableSeq');
	use_ok('Bio::AlignIO');
}

my $DEBUG = test_debug();

my $aa_align = Bio::SimpleAlign->new();
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

is $dna_aln->length, 36;
is $dna_aln->num_residues, 99;
is $dna_aln->num_sequences, 3;
is $dna_aln->consensus_string(50), "atgctgat?gacgtacgc????cgctagcact?aga";

$dna_aln->verbose(-1);
my $replicates;
ok $replicates = &bootstrap_replicates($dna_aln,3);

is scalar @$replicates, 3;
my $repl_aln = pop @$replicates;
is $repl_aln->num_sequences, 3;

##use IO::String;
##use Bio::AlignIO;
##my $string;
##my $out = IO::String->new($string);
##
##my $strout = Bio::AlignIO->new(-fh   => $out,'-format' => 'pfam');
##$strout->write_aln($repl_aln);
##is $string, "";
