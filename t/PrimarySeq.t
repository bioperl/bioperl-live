# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 53);
	
    use_ok('Bio::PrimarySeq');
    use_ok('Bio::Location::Simple');
    use_ok('Bio::Location::Fuzzy');
    use_ok('Bio::Location::Split');
}

my $seq = Bio::PrimarySeq->new(
					 '-seq'            => 'TTGGTGGCGTCAACT',
			       '-display_id'       => 'new-id',
			       '-alphabet'         => 'dna',
			       '-accession_number' => 'X677667',
			       '-desc'             => 'Sample Bio::Seq object');
ok defined $seq;
isa_ok $seq,'Bio::PrimarySeqI';
is $seq->accession_number(), 'X677667';
is $seq->seq(), 'TTGGTGGCGTCAACT';
is $seq->display_id(), 'new-id';
is $seq->alphabet(), 'dna';
is $seq->is_circular(), undef;
ok $seq->is_circular(1);
is $seq->is_circular(0), 0;

# check IdentifiableI and DescribableI interfaces
isa_ok $seq,'Bio::IdentifiableI';
isa_ok $seq,'Bio::DescribableI';
# make sure all methods are implemented
is $seq->authority("bioperl.org"), "bioperl.org";
is $seq->namespace("t"), "t";
is $seq->version(0), 0;
is $seq->lsid_string(), "bioperl.org:t:X677667";
is $seq->namespace_string(), "t:X677667.0";
is $seq->description(), 'Sample Bio::Seq object';
is $seq->display_name(), "new-id";

my $location = Bio::Location::Simple->new('-start' => 2, 
													  '-end' => 5,
													  '-strand' => -1);
is ($seq->subseq($location), 'ACCA');

my $splitlocation = Bio::Location::Split->new();
$splitlocation->add_sub_Location( Bio::Location::Simple->new(
								 '-start' => 1,
							    '-end'   => 4,
							    '-strand' => 1));

$splitlocation->add_sub_Location( Bio::Location::Simple->new(
                         '-start' => 7,
							    '-end'   => 12,
							    '-strand' => -1));

is( $seq->subseq($splitlocation), 'TTGGTGACGC');

my $fuzzy = Bio::Location::Fuzzy->new(-start => '<3',
												 -end   => '8',
												 -strand => 1);

is( $seq->subseq($fuzzy), 'GGTGGC');

my $trunc = $seq->trunc(1,4);
isa_ok $trunc, 'Bio::PrimarySeqI';
is $trunc->seq(), 'TTGG' or diag("Expecting TTGG. Got ".$trunc->seq());

$trunc = $seq->trunc($splitlocation);
isa_ok($trunc, 'Bio::PrimarySeqI');
is( $trunc->seq(), 'TTGGTGACGC');

$trunc = $seq->trunc($fuzzy);
isa_ok($trunc, 'Bio::PrimarySeqI');
is( $trunc->seq(), 'GGTGGC');

my $rev = $seq->revcom();
isa_ok($rev, 'Bio::PrimarySeqI');

is $rev->seq(), 'AGTTGACGCCACCAA' or diag('revcom() failed, was ' . $rev->seq());

#
# Translate
#

my $aa = $seq->translate(); # TTG GTG GCG TCA ACT
is $aa->seq, 'LVAST', "Translation: ". $aa->seq;

# tests for non-standard initiator codon coding for
# M by making translate() look for an initiator codon and
# terminator codon ("complete", the 5th argument below)
$seq->seq('TTGGTGGCGTCAACTTAA'); # TTG GTG GCG TCA ACT TAA
$aa = $seq->translate(undef, undef, undef, undef, 1);
is $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# same test as previous, but using named parameter
$aa = $seq->translate(-complete => 1);
is $aa->seq, 'MVAST', "Translation: ". $aa->seq;

# find ORF, ignore codons outside the ORF or CDS
$seq->seq('TTTTATGGTGGCGTCAACTTAATTT'); # ATG GTG GCG TCA ACT
$aa = $seq->translate(-orf => 1);
is $aa->seq, 'MVAST*', "Translation: ". $aa->seq;

# smallest possible ORF
$seq->seq("ggggggatgtagcccc"); # atg tga
$aa = $seq->translate(-orf => 1);
is $aa->seq, 'M*', "Translation: ". $aa->seq;

# same as previous but complete, so * is removed
$aa = $seq->translate(-orf => 1,
                      -complete => 1);
is $aa->seq, 'M', "Translation: ". $aa->seq;

# ORF without termination codon
# should warn, let's change it into throw for testing
$seq->verbose(2);
$seq->seq("ggggggatgtggcccc"); # atg tgg ccc
eval { $seq->translate(-orf => 1); };
if ($@) {
    like( $@, qr/atgtggcccc\n/);
	$seq->verbose(-1);
	$aa = $seq->translate(-orf => 1);
    is $aa->seq, 'MWP', "Translation: ". $aa->seq;
}
$seq->verbose(0);

# use non-standard codon table where terminator is read as Q
$seq->seq('ATGGTGGCGTCAACTTAG'); # ATG GTG GCG TCA ACT TAG
$aa = $seq->translate(-codontable_id => 6);
is $aa->seq, 'MVASTQ' or diag("Translation: ". $aa->seq);

# insert an odd character instead of terminating with *
$aa = $seq->translate(-terminator => 'X');
is $aa->seq, 'MVASTX' or diag("Translation: ". $aa->seq);

# change frame from default
$aa = $seq->translate(-frame => 1); # TGG TGG CGT CAA CTT AG
is $aa->seq, 'WWRQL' or diag("Translation: ". $aa->seq);

$aa = $seq->translate(-frame => 2); # GGT GGC GTC AAC TTA G
is $aa->seq, 'GGVNL' or diag("Translation: ". $aa->seq);

# TTG is initiator in Standard codon table? Afraid so.
$seq->seq("ggggggttgtagcccc"); # ttg tag
$aa = $seq->translate(-orf => 1);
is $aa->seq, 'L*' or diag("Translation: ". $aa->seq);

# Replace L at 1st position with M by setting complete to 1 
$seq->seq("ggggggttgtagcccc"); # ttg tag
$aa = $seq->translate(-orf => 1,
							 -complete => 1);
is $aa->seq, 'M' or diag("Translation: ". $aa->seq);

# Ignore non-ATG initiators (e.g. TTG) in codon table
$seq->seq("ggggggttgatgtagcccc"); # atg tag
$aa = $seq->translate(-orf => 1,
							 -start => "atg",
							 -complete => 1);
is $aa->seq, 'M' or diag("Translation: ". $aa->seq);


# test for character '?' in the sequence string
is $seq->seq('TTGGTGGCG?CAACT'), 'TTGGTGGCG?CAACT';

# test for some aliases
$seq = Bio::PrimarySeq->new(-id          => 'aliasid',
									 -description => 'Alias desc');
is($seq->description, 'Alias desc');
is($seq->display_id, 'aliasid');

# test that x's are ignored and n's are assumed to be 'dna'
$seq->seq('atgxxxxxx');
is($seq->alphabet,'dna');
$seq->seq('atgnnnnnn');
is($seq->alphabet,'dna');
