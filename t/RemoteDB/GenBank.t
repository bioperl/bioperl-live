# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;

	test_begin(-tests => 44,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);

	use_ok('Bio::DB::GenBank');
}

my %expected_lengths = (
    'MUSIGHBA1' => 408,
    'AF303112'  => 1611,
    'AF303112.1' => 1611,
    'AF041456'  => 1156,
    'CELRABGDI' => 1743,
    'JH374761'  => 38055
);

my ($gb, $seq, $seqio, $seqin);

#
# Bio::DB::GenBank
#

ok $gb = Bio::DB::GenBank->new('-delay'=>0), 'Bio::DB::GenBank';

# get a single seq
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('J00522');1};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Do you have network access? Skipping GenBank tests: $@", 4 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests: $@", 3 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_version('AF303112.1');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests: $@", 2 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_gi('405830');};
    skip "Couldn't connect to Genbank with Bio::DB::GenBank.pm. Transient network problems? Skipping GenBank tests: $@", 1 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
}

$seq = $seqio = undef;

# batch mode
SKIP: {
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Batch access test failed for Genbank. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id}, $s->display_id;
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test the temporary file creation and fasta
ok $gb = Bio::DB::GenBank->new('-format' => 'fasta', '-retrievaltype' => 'tempfile', '-delay' => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('J00522');};
    skip "Couldn't connect to complete GenBank tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 6 if $@;
    # last part of id holds the key
    is $seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to complete GenBank tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 5 if $@;
    # last part of id holds the key
    is $seq->length, $expected_lengths{(split(/\|/,$seq->display_id))[-1]}, $seq->display_id;
    # batch mode requires genbank format
    $gb->request_format("gb");
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Couldn't connect to complete GenBank batch tests with a tempfile with Bio::DB::GenBank.pm. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id};
        undef $gb; # test the case where the db is gone,
        # but a temp file should remain until seqio goes away.
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test pipeline creation
ok $gb = Bio::DB::GenBank->new('-retrievaltype' => 'pipeline', '-delay' => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('J00522');};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 6 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('AF303112');};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 5 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seqio = $gb->get_Stream_by_id([qw(J00522 AF303112 2981014)]);};
    skip "Couldn't connect to complete GenBank tests with a pipeline with Bio::DB::GenBank.pm. Skipping those tests", 4 if $@;
    my $done = 0;
    while (my $s = $seqio->next_seq) {
        is $s->length, $expected_lengths{$s->display_id}, $s->display_id;
        undef $gb; # test the case where the db is gone,
        # but the pipeline should remain until seqio goes away
        $done++;
    }
    skip('No seqs returned', 4) if !$done;
    is $done, 3;
}

$seq = $seqio = undef;

# test contig retrieval
ok $gb = Bio::DB::GenBank->new('-delay'  => 0, '-format' => 'gbwithparts');
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('JH374761');};
    skip "Couldn't connect to GenBank with Bio::DB::GenBank.pm. Skipping those tests", 3 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    # now to check that postprocess_data in NCBIHelper catches CONTIG...
    ok $gb = Bio::DB::GenBank->new('-delay' => 0, '-format' => 'gb');
    eval {$seq = $gb->get_Seq_by_id('JH374761');};
    skip "Couldn't connect to GenBank with Bio::DB::GenBank.pm. Skipping those tests", 1 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
}

$seq = $seqio = undef;

# bug 1405
my @result;
ok $gb = Bio::DB::GenBank->new(-format => 'Fasta', -seq_start  => 2, -seq_stop   => 7);
SKIP: {
    eval {$seq = $gb->get_Seq_by_acc("A11111");};
    skip "Couldn't connect to complete GenBank tests. Skipping those tests", 15 if $@;
    is $seq->length, 6;
    # complexity tests
    ok $gb = Bio::DB::GenBank->new(-format => 'Fasta', -complexity => 0);
    eval {$seqin = $gb->get_Stream_by_acc("5");};
    skip "Couldn't connect to complete GenBank tests. Skipping those tests", 13 if $@;
    @result = (1136, 'dna', 342, 'protein');
    while ($seq = $seqin->next_seq) {
        is $seq->length, shift(@result);
        is $seq->alphabet, shift(@result);
    }
    is @result, 0;
    # Real batch retrieval using epost/efetch
    # these tests may change if integrated further into Bio::DB::Gen*
    # Currently only useful for retrieving GI's via get_seq_stream
    $gb = Bio::DB::GenBank->new();
    eval {$seqin = $gb->get_seq_stream(-uids => [4887706 ,431229, 147460], -mode => 'batch');};
    skip "Couldn't connect to complete GenBank batchmode epost/efetch tests. Skipping those tests", 8 if $@;
    my %result = ('M59757' => 12611 ,'X76083'=> 3140, 'J01670'=> 1593);
	my $ct = 0;
    while ($seq = $seqin->next_seq) {
		$ct++;
		my $acc = $seq->accession;
        ok exists $result{ $acc };
        is $seq->length, $result{ $acc };
		delete $result{$acc};
    }
    skip('No seqs returned', 8) if !$ct;
	is $ct, 3;
    is %result, 0;
}
