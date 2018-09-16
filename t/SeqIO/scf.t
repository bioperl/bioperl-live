# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 80);
	
    use_ok('Bio::SeqIO::scf');
    use_ok('Bio::Seq::SequenceTrace');
}

my $verbose = test_debug();

ok my $in_scf = Bio::SeqIO->new(-file => test_input_file('chad100.scf'),
				-format => 'scf',
				-verbose => $verbose);

my $swq = $in_scf->next_seq();

isa_ok($swq,"Bio::Seq::SequenceTrace");

cmp_ok (length($swq->seq()), '>', 10);
my $qualities = join(' ',@{$swq->qual()});

cmp_ok (length($qualities), '>', 10);
my $id = $swq->id();
is ($swq->id(), "ML4942R");

my $a_channel = $swq->trace("a");
cmp_ok (scalar(@$a_channel), '>', 10);
my $c_channel = $swq->trace("c");
cmp_ok (scalar(@$c_channel), '>', 10);
my $g_channel = $swq->trace("g");
cmp_ok (scalar(@$g_channel), '>', 10);
my $t_channel = $swq->trace("t");
cmp_ok (scalar(@$t_channel), '>', 10);

my $ref = $swq->peak_indices();
my @indices = @$ref;
my $indexcount = 761;
is (scalar(@indices), $indexcount);

#use Data::Dumper;
#----------------------------------------
isa_ok $swq->seq_obj, 'Bio::Seq::Quality';
isa_ok $swq->qual_obj, 'Bio::Seq::Quality';
is $swq->alphabet, 'dna', 'alphabet';

is $swq->display_id, 'ML4942R', 'display_id';
like $swq->primary_id, qr/HASH/, 'primary_id is the stringified memory position';
is $swq->primary_id('ABC'), 'ABC', 'set primary_id';

is $swq->accession_number, 'unknown', 'accession_number';
is $swq->desc, undef, 'desc';
is $swq->desc('test'), 'test', 'desc';
is $swq->id, 'ML4942R', 'id';
is $swq->id('test'), 'test', 'id';
is length($swq->seq), $indexcount, 'seq';



my $len = 7; 
my $start = $swq->length-$len+1;
my $end = $swq->length;

is $swq->subseq($start,$end), 'cctcaag', 'subseq';
is $swq->baseat($start), 'c', 'baseat';
is $swq->qualat($start), '18', 'qualat';

is $swq->trace_value_at('a',$start), '482', 'trace_value_at';

TODO: {
    local $TODO = 'documentation and code for accuracies() do not match' if 1;
    is $swq->accuracies('a',$start), '482', 'accuracies';
}
my $qualstring = join(' ',@{$swq->subqual($start,$end)});
is ($qualstring, '18 18 21 15 8 8 8');

my $refs = $swq->sub_peak_index($start,$end);
is @$refs, $len, 'sub_peak_index';
is $swq->peak_index_at($start), 8819, 'peak_index_at';

my $indices_at_end = join(' ',@{$swq->sub_peak_index($start,$end)});
is($indices_at_end, '8819 8831 8843 8853 8862 8873 8891');

my $swq_end = $swq->trace_length();
my $swq_start = $swq_end - $len +1;
my $subtrace_a = join(' ',@{$swq->sub_trace('a',$swq_start,$swq_end)});
is $subtrace_a, '13 3 0 0 75 274 431';

my $swq2 = $swq->sub_trace_object(1,5);
#$traces2->verbose(-1);

isa_ok($swq2, 'Bio::Seq::SequenceTrace');

$swq2->_synthesize_traces(), 1; # this should not be a private method! Heikki
$swq2->set_accuracies(), 1;

is $swq->accuracy_at('a',1), '755', 'accuracy_at';

#----------------------------------------


warn("Now checking version3...\n") if $verbose;
my $in_scf_v3 = Bio::SeqIO->new(-file => test_input_file('version3.scf'),
				-format => 'scf',
				-verbose => $verbose);

my $v3 = $in_scf_v3->next_seq();
isa_ok($v3, 'Bio::Seq::SequenceTrace');
my $ind = $v3->peak_indices();
my @ff = @$ind;

@indices = @{$v3->peak_indices()};
is (scalar(@indices), 1106);

my %header = %{$in_scf_v3->get_header()};
is $header{bases}, 1106;
is $header{samples},  14107;

# is the Bio::Seq::SequenceTrace AnnotatableI?
my $ac = $v3->annotation();

isa_ok($ac,"Bio::Annotation::Collection");

my @name_comments = grep {$_->tagname() eq 'NAME'} 
$ac->get_Annotations('comment');

is $name_comments[0]->as_text(), 'Comment: IIABP1D4373';

# also get comments this way...
$ac = $in_scf_v3->get_comments();

isa_ok($ac,"Bio::Annotation::Collection");

@name_comments = grep {$_->tagname() eq 'NAME'} 
$ac->get_Annotations('comment');

is $name_comments[0]->as_text(), 'Comment: IIABP1D4373';

my @conv_comments = grep {$_->tagname() eq 'CONV'} 
$ac->get_Annotations('comment');

is $conv_comments[0]->as_text(), 'Comment: phred version=0.990722.h';

# is the SequenceTrace object annotated?
my $st_ac = $swq->annotation();

isa_ok ($st_ac, "Bio::Annotation::Collection");

my @ann =   $st_ac->get_Annotations();

is $ann[0]->tagname, 'SIGN';
is $ann[2]->text, 'SRC3700';
is $ann[5]->tagname, 'LANE';
is $ann[5]->text, 89;
is $ann[6]->text, 'phred version=0.980904.e';
is $ann[8]->text, 'ABI 373A or 377';

my $outfile = test_output_file();
my $out_scf = Bio::SeqIO->new(-file => ">$outfile",
			      -format => 'scf',
			      -verbose => $verbose);

# Bug 2196 - commentless scf

my $in = Bio::SeqIO->new(-file => test_input_file('13-pilE-F.scf'),
			 -format => 'scf',
			 -verbose => $verbose);

my $seq = $in->next_seq;

ok ($seq);

isa_ok($seq, 'Bio::Seq::SequenceTrace');

$ac = $seq->annotation;

isa_ok($ac, 'Bio::Annotation::Collection');

@name_comments = grep {$_->tagname() eq 'NAME'} 
$ac->get_Annotations('comment');

is $name_comments[0], undef;

@conv_comments = grep {$_->tagname() eq 'CONV'} 
$ac->get_Annotations('comment');

is $conv_comments[0], undef;

# the new way

warn("Now testing the _writing_ of scfs\n") if $verbose;

$out_scf->write_seq(-target	=>	$v3,
		    -MACH		=>	'CSM sequence-o-matic 5000',
		    -TPSW		=>	'trace processing software',
		    -BCSW		=>	'basecalling software',
		    -DATF		=>	'AM_Version=2.00',
		    -DATN		=>	'a22c.alf',
		    -CONV		=>	'Bioperl-scf.pm');

ok( -s $outfile && ! -z "$outfile" );

# TODO? tests below don't do much

$out_scf = Bio::SeqIO->new(-verbose => 1,
			   -file => ">$outfile",
			   -format => 'scf');

$swq = Bio::Seq::Quality->new(-seq =>'ATCGATCGAA',
			      -qual =>"10 20 30 40 50 20 10 30 40 50",
			      -alphabet =>'dna');

my $trace = Bio::Seq::SequenceTrace->new(-swq => $swq);

$out_scf->write_seq(	-target	        =>	$trace,
			-MACH		=>	'CSM sequence-o-matic 5000',
			-TPSW		=>	'trace processing software',
			-BCSW		=>	'basecalling software',
			-DATF		=>	'AM_Version=2.00',
			-DATN		=>	'a22c.alf',
			-CONV		=>	'Bioperl-scf.pm' );

warn("Trying to write an scf with a subset of a real scf...\n") if $verbose;
$out_scf = Bio::SeqIO->new(-verbose => 1,
			   -file => ">$outfile",
			   -format => 'scf');

$in_scf_v3 = Bio::SeqIO->new(-file => test_input_file('version3.scf'),
			     -format => 'scf',
			     -verbose => $verbose);
$v3 = $in_scf_v3->next_seq();

my $sub_v3 = $v3->sub_trace_object(5,50);

#warn("The subtrace object is this:\n") if $DEBUG;

$out_scf->write_seq(-target => $sub_v3 );

my $in_scf_v2 = Bio::SeqIO->new(-file => test_input_file('version2.scf'),
				-format => 'scf',
				-verbose => $verbose);
$v3 = $in_scf_v2->next_seq();
ok($v3);

$out_scf = Bio::SeqIO->new(-file   => ">$outfile",
                           -format => "scf");
$out_scf->write_seq( -target  => $v3,
                     -version => 2 );

# simple round trip tests (bug 2881)

my %file_map = (
	# filename         # write_seq args
	'chad100.scf'	=> 1,
	'13-pilE-F.scf' => 1,
	'version2.scf'  => 1,
	'version3.scf'  => 1
	);

for my $f (sort keys %file_map) {
	my $outfile = test_output_file();
	my $in = Bio::SeqIO->new(-file => test_input_file($f),
							 -format => 'scf');
	my $out = Bio::SeqIO->new(-file => ">$outfile",
							 -format => 'scf');
	
	my $seq1 = $in->next_seq();
	isa_ok($seq1, 'Bio::Seq::SequenceTrace');
	
	ok($out->write_seq(-target => $seq1));
	
	my $in2 = Bio::SeqIO->new(-file => "<$outfile",
							  -format => 'scf');
	my $seq2 = $in2->next_seq();
	isa_ok($seq2, 'Bio::Seq::SequenceTrace');
	if ($seq1->display_id) {
		TODO: {
			local $TODO = "display_id doesn't round trip yet";
			is($seq1->display_id, $seq2->display_id, 'display_id matches');
		}
	}
	is_deeply($seq1->qual, $seq2->qual, 'qual scores match');
}

# synthesizing traces roundtrip (bug 2881):

my @sequences=('A',
               'ATGGAGCTCATCAAAGAATCGACTCATATATCCATCCCTGAACGGCTGACTCACATTAATGGTTGA');
foreach my $sequence (@sequences) {
    my $qualstr=join ' ', map { 65 } (1 .. length($sequence));
    my $seq_qual=Bio::Seq::Quality->new(-seq=>$sequence, -qual=>$qualstr);

    my $outfile=test_output_file();
    my $out=Bio::SeqIO->new(-file=>">$outfile", -format=>'scf');
    $out->write_seq(-target=>$seq_qual);

    my $in=Bio::SeqIO->new(-file=>$outfile, -format=>'scf');
    my $in_seq=$in->next_seq();

    is_deeply($seq_qual, $in_seq->{swq}, 'Bio::Sequence::Quality matches');
}

