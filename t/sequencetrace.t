# -*-Perl-*-

use strict;
use Dumpvalue;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;
use Bio::Seq::Quality;

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 5;
}

my $dumper = new Dumpvalue();
$dumper->veryCompact(1);
# $dumper->compactDump(1);
my $DEBUG = $ENV{'BIOPERLDEBUG'};

        # redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");

# print("Checking if the Bio::Seq::SequenceTrace module could be used...\n") if $DEBUG;
        # test 1
use Bio::Seq::SequenceTrace;
     # test 1
ok(1);

# print("Reading an scf...\n");

my $in_scf_v3 = Bio::SeqIO->new('-file' => Bio::Root::IO->catfile
                    ("t","data",
                     "version3.scf"),
                    '-format' => 'scf',
                    );

my $trace = $in_scf_v3->next_seq();

# print("Testing those values...\n");
# print("Length: ".$trace->length()."\n");
     # at the very end
my $start = $trace->length()-19;
my $end = $trace->length();

# print("Testing subseq from the end...".$trace->subseq($start,$end)."\n");
     # test 2
ok ($trace->subseq($trace->length()-19,$trace->length()) eq "CCCCTTTCCCAACAGCACCG");
# print("Testing the qualities for those bases...".join(' ',@{$trace->subqual($start,$end)})."\n");
my $qualstring = join(' ',@{$trace->subqual($start,$end)});
     # test 3
ok ($qualstring eq "12 10 7 7 9 7 7 9 13 9 9 9 6 6 6 8 8 8 6 6");
# print("Testing getting sub qual indices\n");

# print("Testing the trace indices for those bases ($start->$end):...\n");
my $ref = $trace->sub_peak_index($start,$end);
my @temp = @{$ref};

my $indices_at_end = join(' ',@{$trace->sub_peak_index($start,$end)});
     # test 4
ok($indices_at_end eq "13863 13874 13883 13898 13905 13922 13934 13952 13966 13975 13982 14003 14013 14026 14037 14056 14061 14084 14093 14099");
# print("Getting all of the trace values for that range\n");
my $trace_end = $trace->trace_length();
my $trace_start = $trace_end - 19;
my $subtrace_a = join(' ',@{$trace->sub_trace('a',$trace_start,$trace_end)});
     # test 5
(ok $subtrace_a eq "63 61 68 82 101 120 135 145 148 143 131 111 85 59 37 18 4 0 3 6");
     # print("scf_dump ing...\n");
     # $trace->scf_dump();
     # print("The traces are:\n");
     # $trace->_dump_traces();



# whew! now given a subset of bases, get their traces....
my $traces2 = $trace->sub_trace_object(1,5);
$traces2->verbose(-1);


# print("Attempting to synthesize traces for this object:\n");
# print("The sequence is : ".$traces2->seq()."\n");
# print("The qualities are : ".join(' ',@{$traces2->qual()})."\n");
# print("The length is : ".$traces2->length()."\n");
# $dumper->dumpValue($traces2);
# can you synthesize false traces?

$traces2->_synthesize_traces();


$traces2->set_accuracies();
print("This is an scf dump:\n");
$traces2->scf_dump();

sub the_old_way {
        my $start2 = 1;
        my $end2 = 5;
        my $subtraces;
        $subtraces->{base_start} = $start2;
        $subtraces->{base_end} = $end2;
        $subtraces->{sequence} = $trace->subseq($start2,$end2);
        $subtraces->{qualities} = join(' ',$trace->subqual($start2,$end2));
        $subtraces->{indices} = $trace->sub_trace_index($start2,$end2);
        my @temp = @{$subtraces->{indices}};
        $subtraces->{trace_start} = $temp[0];
        $subtraces->{trace_end} = $temp[$#temp];
        foreach (qw(a t g c)) {
             $subtraces->{traces}->{$_} = $trace->sub_trace($_,$subtraces->{trace_start},$subtraces->{trace_end});
        }

        $dumper->dumpValue($subtraces);
}


