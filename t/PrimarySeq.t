# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
    plan tests => 33;
}
use Bio::PrimarySeq;
use Bio::Location::Simple;
use Bio::Location::Fuzzy;
use Bio::Location::Split;

ok(1);

my $seq = Bio::PrimarySeq->new('-seq'              =>'TTGGTGGCGTCAACT',
			       '-display_id'       => 'new-id',
			       '-alphabet'          => 'dna',
			       '-accession_number' => 'X677667',
			       '-desc'             =>'Sample Bio::Seq object');
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

my $location = new Bio::Location::Simple('-start' => 2, '-end' => 5,
					 '-strand' => -1);
ok ($seq->subseq($location), 'ACCA');

my $splitlocation = new Bio::Location::Split();
$splitlocation->add_sub_Location( new Bio::Location::Simple('-start' => 1,
							    '-end'   => 4,
							    '-strand' => 1));

$splitlocation->add_sub_Location( new Bio::Location::Simple('-start' => 7,
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

ok $rev->seq(), 'AGTTGACGCCACCAA', 'rev faile was ' . $rev->seq();

#
# Translate
#

my $aa = $seq->translate();
ok $aa->seq, 'LVAST', "Translation: ". $aa->seq;

$seq->seq('TTGGTGGCGTCAACTTAA');
$aa = $seq->translate(undef, undef, undef, undef, 1);

# tests for non-Methionin initiator codon (AGT) coding for M
ok $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# test for character '?' in the sequence string
ok $seq->seq('TTGGTGGCG?CAACT'), 'TTGGTGGCG?CAACT';
