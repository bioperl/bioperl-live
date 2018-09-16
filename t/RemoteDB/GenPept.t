# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
	use Bio::Root::Test;

	test_begin(-tests => 21,
			   -requires_modules => [qw(IO::String
									    LWP::UserAgent
										HTTP::Request::Common)],
			   -requires_networking => 1);

	use_ok('Bio::DB::GenPept');
}

my %expected_lengths = (
    'AAC06201'  => 353,
    'AAD15290'  => 136,
    'P31383'    => 635,
    '2AAA_YEAST' => 635
);

my ($gb, $seq, $seqio);

#
# Bio::DB::GenPept
#
ok $gb = Bio::DB::GenPept->new();
SKIP: {
    eval {$seqio = $gb->get_seq_stream(-uids => [2981015, 1621261, 195055], -mode => 'batch');};
    skip "Couldn't connect to complete GenPept tests. Skipping those tests", 8 if $@;
    my %result = ('AAC06201' => 353, 'CAB02640' => 193, 'AAD15290' => 136);
    my $ct = 0;
    while ($seq = $seqio->next_seq) {
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

$seq = $seqio = undef;

ok $gb = Bio::DB::GenPept->new('-delay' => 0);
SKIP: {
    eval {$seq = $gb->get_Seq_by_id('195055');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 10 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seq = $gb->get_Seq_by_acc('AAC06201');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 9 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;
    eval {$seqio = $gb->get_Stream_by_id([qw(AAC06201 195055)]);};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 8 if $@;
    my $done = 0;
    while( my $s = $seqio->next_seq ) {
        is $s->length, $expected_lengths{$s->display_id}, $s->display_id;
        $done++;
    }
    skip('No seqs returned', 8) if !$done;
    is $done, 2;
    # swissprot genpept parsing
    eval {$seq = $gb->get_Seq_by_acc('P31383');};
    skip "Couldn't connect to Genbank with Bio::DB::GenPept.pm. Skipping those tests", 5 if $@;
    is $seq->length, $expected_lengths{$seq->display_id}, $seq->display_id;

    # test dbsource stuff
    # small chance this might change but hopefully not
    my @annot = $seq->annotation->get_Annotations('dblink');
    cmp_ok(scalar(@annot), '>', 31);
    is $annot[0]->database, 'UniProtKB';
    is $annot[0]->primary_id, '2AAA_YEAST';
    is (($seq->annotation->get_Annotations('swissprot_dates'))[0]->value, 'Jul 1, 1993');
}
