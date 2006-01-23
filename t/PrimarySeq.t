# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl t/test.t'

use strict;
BEGIN { 
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
	eval { require Test; };
	if( $@ ) {
		use lib 't';
	}
	use Test;
	plan tests => 48;
}
use Bio::PrimarySeq;
use Bio::Location::Simple;
use Bio::Location::Fuzzy;
use Bio::Location::Split;

ok(1);

my $seq = Bio::PrimarySeq->new(
					 '-seq'              => 'TTGGTGGCGTCAACT',
			       '-display_id'       => 'new-id',
			       '-alphabet'         => 'dna',
			       '-accession_number' => 'X677667',
			       '-desc'             => 'Sample Bio::Seq object');
ok defined $seq;
ok $seq->isa('Bio::PrimarySeqI');
ok $seq->accession_number(), 'X677667';
ok $seq->seq(), 'TTGGTGGCGTCAACT';
ok $seq->display_id(), 'new-id';
ok $seq->alphabet(), 'dna';
ok $seq->is_circular(), undef;
ok $seq->is_circular(1);
ok $seq->is_circular(0), 0;

# check IdentifiableI and DescribableI interfaces
ok $seq->isa('Bio::IdentifiableI');
ok $seq->isa('Bio::DescribableI');
# make sure all methods are implemented
ok $seq->authority("bioperl.org"), "bioperl.org";
ok $seq->namespace("t"), "t";
ok $seq->version(0), 0;
ok $seq->lsid_string(), "bioperl.org:t:X677667";
ok $seq->namespace_string(), "t:X677667.0";
ok $seq->description(), 'Sample Bio::Seq object';
ok $seq->display_name(), "new-id";

my $location = new Bio::Location::Simple('-start' => 2, 
													  '-end' => 5,
													  '-strand' => -1);
ok ($seq->subseq($location), 'ACCA');

my $splitlocation = new Bio::Location::Split();
$splitlocation->add_sub_Location( new Bio::Location::Simple(
								 '-start' => 1,
							    '-end'   => 4,
							    '-strand' => 1));

$splitlocation->add_sub_Location( new Bio::Location::Simple(
                         '-start' => 7,
							    '-end'   => 12,
							    '-strand' => -1));

ok( $seq->subseq($splitlocation), 'TTGGTGACGC');

my $fuzzy = new Bio::Location::Fuzzy(-start => '<3',
												 -end   => '8',
												 -strand => 1);

ok( $seq->subseq($fuzzy), 'GGTGGC');

my $trunc = $seq->trunc(1,4);
ok defined $trunc;
ok $trunc->seq(), 'TTGG', "Expecting TTGG. Got ".$trunc->seq();

$trunc = $seq->trunc($splitlocation);
ok( defined $trunc);
ok( $trunc->seq(), 'TTGGTGACGC');

$trunc = $seq->trunc($fuzzy);
ok( defined $trunc);
ok( $trunc->seq(), 'GGTGGC');

my $rev = $seq->revcom();
ok defined $rev; 

ok $rev->seq(), 'AGTTGACGCCACCAA', 'revcom() failed, was ' . $rev->seq();

#
# Translate
#

my $aa = $seq->translate(); # TTG GTG GCG TCA ACT
ok $aa->seq, 'LVAST', "Translation: ". $aa->seq;

# tests for non-standard initiator codon coding for
# M by making translate() look for an initiator codon and
# terminator codon ("complete", the 5th argument below)
$seq->seq('TTGGTGGCGTCAACTTAA'); # TTG GTG GCG TCA ACT TAA
$aa = $seq->translate(undef, undef, undef, undef, 1);
ok $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# same test as previous, but using named parameter
$aa = $seq->translate(-complete => 1);
ok $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# find ORF, ignore codons outside the ORF or CDS
$seq->seq('TTTTATGGTGGCGTCAACTTAATTT'); # ATG GTG GCG TCA ACT
$aa = $seq->translate(-orf => 1);
ok $aa->seq, 'MVAST*', "Translation: ". $aa->seq;

# smallest possible ORF
$seq->seq("ggggggatgtagcccc"); # atg tga
$aa = $seq->translate(-orf => 1);
ok $aa->seq, 'M*', "Translation: ". $aa->seq;

# same as previous but complete, so * is removed
$aa = $seq->translate(-orf => 1,
                      -complete => 1);
ok $aa->seq, 'M', "Translation: ". $aa->seq;

# ORF without termination codon
# should warn, let's change it into throw for testing
$seq->verbose(2);
$seq->seq("ggggggatgtggcccc"); # atg tgg ccc

eval {$aa = $seq->translate(-orf => 1);};
if ($@) {
    ok  1 if $@ =~ /atgtggcccc\n/;
    #ok $aa->seq, 'MWP', "Translation: ". $aa->seq; 
}
$seq->verbose(0);

# use non-standard codon table where terminator is read as Q
$seq->seq('ATGGTGGCGTCAACTTAG'); # ATG GTG GCG TCA ACT TAG
$aa = $seq->translate(-codontable_id => 6);
ok $aa->seq, 'MVASTQ', "Translation: ". $aa->seq;

# insert an odd character instead of terminating with *
$aa = $seq->translate(-terminator => 'X');
ok $aa->seq, 'MVASTX', "Translation: ". $aa->seq;

# change frame from default
$aa = $seq->translate(-frame => 1); # TGG TGG CGT CAA CTT AG
ok $aa->seq, 'WWRQL', "Translation: ". $aa->seq;

# TTG is initiator in Standard codon table? Afraid so.
$seq->seq("ggggggttgtagcccc"); # ttg tag
$aa = $seq->translate(-orf => 1);
ok $aa->seq, 'L*', "Translation: ". $aa->seq;

# Replace L at 1st position with M by setting complete to 1 
$seq->seq("ggggggttgtagcccc"); # ttg tag
$aa = $seq->translate(-orf => 1,
							 -complete => 1);
ok $aa->seq, 'M', "Translation: ". $aa->seq;

# Ignore non-ATG initiators (e.g. TTG) in codon table
$seq->seq("ggggggttgatgtagcccc"); # atg tag
$aa = $seq->translate(-orf => 1,
							 -start => "atg",
							 -complete => 1);
ok $aa->seq, 'M', "Translation: ". $aa->seq;



# test for character '?' in the sequence string
ok $seq->seq('TTGGTGGCG?CAACT'), 'TTGGTGGCG?CAACT';

# test for some aliases
$seq = Bio::PrimarySeq->new(-id          => 'aliasid',
									 -description => 'Alias desc');
ok($seq->description, 'Alias desc');
ok($seq->display_id, 'aliasid');

# test that x's are ignored and n's are assumed to be 'dna'
$seq->seq('atgxxxxxx');
ok($seq->alphabet,'dna');
$seq->seq('atgnnnnnn');
ok($seq->alphabet,'dna');
