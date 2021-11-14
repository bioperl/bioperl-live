#!/usr/bin/env perl
use rlib '.';
use strict; use warnings;
use Test::More;
use Helper;

my %notes = (
    'anonymize'   => 'anonymize sequence IDs',
    'composition' => 'base composition',
    'length'      => 'protein sequence length',
    'linearize'   => 'linearize fast sequence',
    'no-gaps'     => 'remove gaps',
    'num-seq'     => 'number of sequences',
    'remove-stop' => 'remove stop codons',
    'revcom'      => 'reverse compliment sequence',
);

# Filenames like
# VS116:7:310... created on MSWindows using --break are invalid
#
if ($^O ne 'MSWin32') {
    $notes{'break'} = 'break into single-sequence files';
}

test_no_arg_opts('bioseq', 'test-bioseq.nuc', \%notes);

my $opts = [
    ['delete', 'order:2', 'delete by order'],
    ['pick', 'order:2', 'pick 1 sequence by order'],
    ['subseq', '10,20', 'get subsequences'],
    ['translate', '1', 'translate DNA'],
    ['reloop', '3', 'reloop a sequence'],
    ];

test_one_arg_opts('bioseq', 'test-bioseq.nuc', $opts);

my $multi_opts = [
    ["--pick order:2,4", 'test-bioseq.nuc',
     'pick-order-2,4.right', 'pick seqs by order delimited by commas'],
    ["--pick order:2-4", 'test-bioseq.nuc',
     'pick-order-2-4.right', 'pick seqs by order with range'],
    ["--restrict EcoRI", 'test-bioseq-re.fas',
     'restrict.right', 'restriction cut'],
    ["--hydroB", 'test-bioseq.pep',
     'hydroB.right', 'Hydrophobicity score'],
    ["--input genbank --output fasta", "test-bioseq.gb",
     "genbank2fast.right", "Genbank => fasta"],
    ["--input genbank --feat2fas", "test-bioseq.gb",
     "feat2fasta.right", "extract gens from a genbank file"]
    ];

    for my $tuple (@$multi_opts) {
	my ($opts, $datafile, $check, $note) = @$tuple;
	run_bio_program('bioseq', $datafile, $opts, $check, {note => $note});
    }


if ($ENV{'BPWRAPPER_INTERNET_TESTS'}) {
    my $multi_opts = [
	["--fetch 'X83553' --output genbank",
	 'X83553.right', 'fetch Genbank file X8553']
	];

    for my $tuple (@$multi_opts) {
	my ($opts, $check, $note) = @$tuple;
	run_bio_program_nocheck('bioseq', '/dev/null', $opts, {note => $note});
    }
}


done_testing();
