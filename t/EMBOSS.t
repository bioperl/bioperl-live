
# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
BEGIN {
    use vars qw($NTESTS);
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $NTESTS = 30;
    plan tests => $NTESTS }

use Bio::Factory::EMBOSS;
use Bio::Root::IO;
use Bio::SeqIO;
use Bio::AlignIO;
my $compseqoutfile = 'dna1.4.compseq';
my $wateroutfile   = 'cysprot.water';
my $consoutfile    = 'cysprot.cons';
END {

    foreach ( $Test::ntest..$NTESTS ) {
	skip("EMBOSS not installed locally or XML::Twig not installed",1);
    }
    unlink($compseqoutfile);
    unlink($wateroutfile);
    unlink($consoutfile);
}

my $verbose = $ENV{'BIOPERLDEBUG'} || -1;
ok(1);

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run.

my $factory = new Bio::Factory::EMBOSS(-verbose => $verbose);

ok($factory);

my $compseqapp = $factory->program('compseq');
if( ! $compseqapp ) { 
    # no EMBOSS installed
    exit();
}

ok($compseqapp);
my %input = ( '-word' => 4,
	      '-sequence' => Bio::Root::IO->catfile('t',
						   'data',
						   'dna1.fa'),
	      '-outfile' => $compseqoutfile);
$compseqapp->run(\%input);
ok(-e $compseqoutfile);
#if( open(IN, $compseqoutfile) ) {
#    while(<IN>) { print }
#}

my $water = $factory->program('water');

ok ($water);

# testing in-memory use of 
my $in = new Bio::SeqIO(-format => 'fasta',
			-file =>  Bio::Root::IO->catfile('t',
							 'data',
							 'cysprot1a.fa'));
my $seq = $in->next_seq();
ok($seq);
my @amino;
$in = new Bio::SeqIO(-format => 'fasta',
		     -file =>  Bio::Root::IO->catfile('t',
						      'data',
						      'amino.fa'));
while( my $s = $in->next_seq) {
    push @amino, $s;
}
$water->run({ '-sequencea' => $seq,
	      '-seqall'    => \@amino,
	      '-gapopen'   => '10.0',
	      '-gapextend' => '0.5',
	      '-outfile'   => $wateroutfile});

ok(-e $wateroutfile);

my $alnin = new Bio::AlignIO(-format => 'emboss',
			    -file   => $wateroutfile);

ok( $alnin);
my $aln = $alnin->next_aln;
ok($aln);
ok($aln->length, 43);
ok($aln->overall_percentage_identity, 100);
ok($aln->average_percentage_identity, 100);

my ($first) = $aln->each_seq();
ok($first->seq(), 'SCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGE');
$aln = $alnin->next_aln;
ok($aln);
ok($aln->length, 339);
ok(sprintf("%.2f",$aln->overall_percentage_identity), 33.04);
ok(sprintf("%.2f",$aln->average_percentage_identity), 40.58);

my $cons = $factory->program('cons');
$cons->verbose(0);
$in = new Bio::AlignIO(-format => 'msf',
		       -file   => Bio::Root::IO->catfile('t',
							 'data',
							 'cysprot.msf'));
my $aln2 = $in->next_aln;
$cons->run({ '-msf'   => $aln2,
	     '-outseq'=> $consoutfile});

ok(-e $consoutfile);


# testing acd parsing and EMBOSSacd methods

$compseqapp = $factory->program('compseq');
exit unless $compseqapp->acd;

ok my $acd = $compseqapp->acd;
ok $compseqapp->acd->name, 'compseq';
ok my $compseq_mand_acd = $compseqapp->acd->mandatory;
ok $compseq_mand_acd->mandatory->qualifier('-word');
ok $compseq_mand_acd->mandatory->qualifier('-supper1'), 0;
ok $acd->qualifier('-ppppppp'), 0;
ok $acd->qualifier('-reverse');
ok $acd->category('-reverse'), 'optional';
ok $acd->values('-reverse'), 'Yes/No';
ok $acd->descr('-reverse'), 'Set this to be true if you also wish to also count words in the reverse complement of a nucleic sequence.';
ok $acd->unnamed('-reverse'), 0;
ok $acd->default('-reverse'), 'No';

## comparing input and ACD qualifiers
## commented out because verbose > 1 prints error messages
## that would confuse users running tests
#
#$compseqapp->verbose(1);
#%input = ( '-word' => 4,
#	   '-outfile' => $compseqoutfile);
#eval {
#    my $a = $compseqapp->run(\%input);
#};
#ok 1 if $@; # '-sequence' missing
#
#%input = ( '-word' => 4,
#	   '-incorrect_option' => 'no value',
#	   '-sequence' => Bio::Root::IO->catfile('t',
#						 'data',
#						 'dna1.fa'),
#	   '-outfile' => $compseqoutfile);
#eval {
#    $compseqapp->run(\%input);
#};
#ok 1 if $@; # -incorrect_option is incorrect
