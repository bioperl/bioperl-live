
use strict;
use Dumpvalue;
use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;
use Bio::Seq::SeqWithQuality;

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 1;
}

my $dumper = new Dumpvalue();
my $DEBUG = $ENV{'BIOPERLDEBUG'};

        # redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");

print("Checking if the Bio::Seq::SequenceTrace module could be used...\n") if $DEBUG;
        # test 1
use Bio::Seq::SequenceTrace;
ok(1);

     # create some random quality object with the same number of qualities and
     # the same identifiers

my $seqobj = Bio::PrimarySeq->new( -seq => "ATCGATCGA",
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                            );
my $string_quals = "10 20 30 40 50 40 30 20 10";
my $indices = "5 10 15 20 25 30 35 40 45";
my $qualobj = Bio::Seq::PrimaryQual->new( -qual => $string_quals,
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                              -trace_indices =>   $indices
                            );

my $swq1 = Bio::Seq::SeqWithQuality->new( -seq	=>	$seqobj,
					-qual		=>	$qualobj);


#     print("\t6b) Testing the getter\n");
#          my @retrieved_indices = @{$swq1->trace_indices()};
#          ok($ti_orig eq join(' ',@retrieved_indices));
#     print("\t6c) Testing the setter\n");
#          my $ti_new = "1 2 3 4 5 6 7 8 9 10";
#          $swq1->trace_indices($ti_new);
#          ok ($ti_new eq join(' ' ,@{$swq1->trace_indices()}));
#     print("\t6d) Testing the subtraceindex at the start (border condition)\n");
#          ok ('1 2 3' eq join(' ',@{$swq1->sub_trace_index(1,3)}));
#     print("\t6d) Testing the subtraceindex at the end (border condition)\n");
#          ok ('7 8 9' eq join(' ',@{$swq1->sub_trace_index(7,9)}));
#     print("\t6d) Testing the subtraceindex in the middle\n");
#          ok ('4 5 6' eq join(' ',@{$swq1->sub_trace_index(4,6)}));
#
