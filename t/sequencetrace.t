# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::SeqIO');
}

my $debug = test_debug();

my $in_scf_v3 = Bio::SeqIO->new('-file' => test_input_file('version3.scf'),
                    '-format' => 'scf');

my $trace = $in_scf_v3->next_seq();

my $start = $trace->length()-19;
my $end = $trace->length();

is ($trace->subseq($trace->length()-19,$trace->length()), "CCCCTTTCCCAACAGCACCG");

my $qualstring = join(' ',@{$trace->subqual($start,$end)});
is ($qualstring, "12 10 7 7 9 7 7 9 13 9 9 9 6 6 6 8 8 8 6 6");

my $ref = $trace->sub_peak_index($start,$end);
my @temp = @{$ref};

my $indices_at_end = join(' ',@{$trace->sub_peak_index($start,$end)});
is($indices_at_end, "13863 13874 13883 13898 13905 13922 13934 13952 13966 13975 13982 14003 14013 14026 14037 14056 14061 14084 14093 14099");

my $trace_end = $trace->trace_length();
my $trace_start = $trace_end - 19;
my $subtrace_a = join(' ',@{$trace->sub_trace('a',$trace_start,$trace_end)});
is $subtrace_a, "63 61 68 82 101 120 135 145 148 143 131 111 85 59 37 18 4 0 3 6";


# whew! now given a subset of bases, get their traces....
my $traces2 = $trace->sub_trace_object(1,5);
$traces2->verbose(-1);

$traces2->_synthesize_traces();

$traces2->set_accuracies();

if ($debug) {
	print("This is an scf dump:\n");
	$traces2->scf_dump();
}
