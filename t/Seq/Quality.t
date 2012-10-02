# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 85);

    use_ok('Bio::Seq::Quality');
}

use Bio::SeqIO;

my $DEBUG = test_debug();

# create some random sequence object with no id
my $seqobj_broken = Bio::Seq::Quality->
    new( -seq => "ATCGATCGA",
	 );

my $seqobj;
lives_ok {
    $seqobj = Bio::Seq::Quality->
	new( -seq => "ATCGATCGA",
	     -id  => 'QualityFragment-12',
	     -accession_number => 'X78121',
	     );
};


# create some random quality object with the same number of qualities
# and the same identifiers
my $string_quals = "10 20 30 40 50 40 30 20 10";
my $qualobj;
lives_ok {
    $qualobj = Bio::Seq::Quality->
	new( -qual => $string_quals,
	     -id  => 'QualityFragment-12',
	     -accession_number => 'X78121',
	     );
};

# check to see what happens when you construct the Quality object
ok my $swq1 = Bio::Seq::Quality->
    new( -seq => "ATCGATCGA",
	 -id  => 'QualityFragment-12',
	 -accession_number => 'X78121',
	 -qual	=>	$string_quals);


print("Testing various weird constructors...\n") if $DEBUG;
print("\ta) No ids, Sequence object, no quality...\n") if $DEBUG;

# w for weird
my $wswq1;
lives_ok {
    $wswq1 = Bio::Seq::Quality->
	new( -seq  => "ATCGATCGA",
	     -qual => "");
};
print $@ if $DEBUG;


print("\tb) No ids, no sequence, quality object...\n") if $DEBUG;
	# note that you must provide a alphabet for this one.
$wswq1 = Bio::Seq::Quality->
    new( -seq => "",
	 -qual => $string_quals,
	 -alphabet => 'dna'
	 );

print("\tc) Absolutely nothing. (HAHAHAHA)...\n") if $DEBUG;
lives_ok {
    $wswq1 = Bio::Seq::Quality->new( -seq => "",
				     -qual => "",
				     -alphabet => 'dna'
				     );
};


print("\td) Absolutely nothing but an ID\n") if $DEBUG;
lives_ok {
    $wswq1 = Bio::Seq::Quality->
	new( -seq => "",
	     -qual => "",
	     -alphabet => 'dna',
	     -id => 'an object with no sequence and no quality but with an id'
	     );
};

print("\td) No sequence, no quality, no ID...\n") if $DEBUG;
warnings_like {
    $wswq1 = Bio::Seq::Quality->
	new( -seq  =>	"",
	     -qual =>	"",
	     -verbose => 0);
} qr/not guess alphabet/i;

print("Testing various methods and behaviors...\n") if $DEBUG;

print("1. Testing the seq() method...\n") if $DEBUG;
	print("\t1a) get\n") if $DEBUG;
	my $original_seq = $swq1->seq();
	is ($original_seq, "ATCGATCGA");
	print("\t1b) set\n") if $DEBUG;
	ok ($swq1->seq("AAAAAAAAAAAA"));
	print("\t1c) get (again, to make sure the set was done.)\n") if $DEBUG;
	is($swq1->seq(), "AAAAAAAAAAAA");
	print("\tSetting the sequence back to the original value...\n") if $DEBUG;
	$swq1->seq($original_seq);


print("2. Testing the qual() method...\n") if $DEBUG;
	print("\t2a) get\n") if $DEBUG;
	my @qual = @{$swq1->qual()};
	my $str_qual = join(' ',@qual);
	is $str_qual, "10 20 30 40 50 40 30 20 10";
	print("\t2b) set\n") if $DEBUG;
	ok $swq1->qual("10 10 10 10 10");
	print("\t2c) get (again, to make sure the set was done.)\n") if $DEBUG;
	my @qual2 = @{$swq1->qual()};
	my $str_qual2 = join(' ',@qual2);
	is($str_qual2, "10 10 10 10 10 0 0 0 0"); ###!
	print("\tSetting the quality back to the original value...\n") if $DEBUG;
	$swq1->qual($str_qual);

print("3. Testing the length() method...\n") if $DEBUG;
	print("\t3a) When lengths are equal...\n") if $DEBUG;
	is($swq1->length(), 9);	
	print("\t3b) When lengths are different\n") if $DEBUG;
	$swq1->qual("10 10 10 10 10");
	isnt ($swq1->length(), "DIFFERENT");


print("6. Testing the subqual() method...\n") if $DEBUG;
     my $t_subqual = "10 20 30 40 50 60 70 80 90";
     $swq1->qual($t_subqual);
     print("\t6d) Testing the subqual at the start (border condition)\n") if $DEBUG;
     ok ('10 20 30' eq join(' ',@{$swq1->subqual(1,3)}));
     print("\t6d) Testing the subqual at the end (border condition)\n") if $DEBUG;
     ok ('70 80 90' eq join(' ',@{$swq1->subqual(7,9)}));
     print("\t6d) Testing the subqual in the middle\n") if $DEBUG;
     ok ('40 50 60' eq join(' ',@{$swq1->subqual(4,6)}));

print("7. Testing cases where quality is zero...\n") if $DEBUG;
$swq1 = Bio::Seq::Quality->new(-seq =>  'G',
                               -qual => '0',
			       );
my $swq2 = Bio::Seq::Quality->new(-seq =>  'G',
                                  -qual => '65',
				  );
is $swq1->length, $swq2->length;

$swq1 = Bio::Seq::Quality->new(-seq =>  'GC',
                               -qual => '0 0',
			       );
$swq2 = Bio::Seq::Quality->new(-seq =>  'GT',
                               -qual => '65 0',
			       );
my $swq3 = Bio::Seq::Quality->new(-seq =>  'AG',
                               -qual => '0 60',
			       );
is $swq1->length, $swq2->length;
is $swq1->length, $swq3->length;


#
# end of test inherited from seqwithquality.t 
#
#################################################################
#
# testing new functionality
#

my $qual = '0 1 2 3 4 5 6 7 8 9 11 12 13';
my $trace = '0 5 10 15 20 25 30 35 40 45 50 55 60';

ok my $seq = Bio::Seq::Quality->new
    ( -qual => $qual,
      -trace_indices => $trace,
      -seq =>  'atcgatcgatcgt',
      -id  => 'human_id',
      -accession_number => 'S000012',
      -verbose => $DEBUG >= 0 ? $DEBUG : 0
);


print("2. Testing the trace() method...\n") if $DEBUG;
	print("\t2a) get\n") if $DEBUG;
	my @trace = @{$seq->trace()};
	my $str_trace = join(' ',@trace);
	is $str_trace, $trace;
	print("\t2b) set\n") if $DEBUG;
	ok $seq->trace("10 10 10 10 10");
	print("\t2c) get (again, to make sure the set was done.)\n") if $DEBUG;
	my @trace2 = @{$seq->trace()};
	my $str_trace2 = join(' ',@trace2);
	is($str_trace2, "10 10 10 10 10 0 0 0 0 0 0 0 0"); ###!
	print("\tSetting the trace back to the original value...\n") if $DEBUG;
	$seq->trace($trace);



is_deeply $seq->qual, [split / /, $qual];
is_deeply $seq->trace, [split / /, $trace];
is_deeply $seq->trace_indices, [split / /, $trace]; #deprecated

is $seq->qual_text, $qual;
is $seq->trace_text, $trace;

is join (' ', @{$seq->subqual(2, 3)}), '1 2';
is $seq->subqual_text(2, 3), '1 2';
is join (' ', @{$seq->subqual(2, 3, "9 9")}), '9 9';
is $seq->subqual_text(2, 3, "8 8"), '8 8';

is join (' ', @{$seq->subtrace(2, 3)}), '5 10';
is $seq->subtrace_text(2, 3), '5 10';
is join (' ', @{$seq->subtrace(2, 3, "9 9")}), '9 9';
is $seq->subtrace_text(2, 3, "8 8"), '8 8';


is $seq->trace_index_at(5), 20;
is join(' ', @{$seq->sub_trace_index(5,6)}), "20 25";

is $seq->baseat(2), 't';
is $seq->baseat(3), 'c';
is $seq->baseat(4), 'g';
is $seq->baseat(5), 'a';


#############################################
#
# same tests using Seq::Meta::Array methods follow ...
#

my $meta = '0 1 2 3 4 5 6 7 8 9 11 12';
$trace = '0 5 10 15 20 25 30 35 40 45 50 55';
my @trace_array = qw(0 5 10 15 20 25 30 35 40 45 50 55);

ok $seq = Bio::Seq::Quality->new
    ( -meta => $meta,
      -seq =>  'atcgatcgatcg',
      -id  => 'human_id',
      -accession_number => 'S000012',
      -verbose => $DEBUG >= 0 ? $DEBUG : 0
);

$seq->named_meta('trace', \@trace_array);

is_deeply $seq->meta, [split / /, $meta];
is_deeply $seq->named_meta('trace'), [split / /, $trace];

is $seq->meta_text, $meta;
is $seq->named_meta_text('trace'), $trace;

is join (' ', @{$seq->submeta(2, 3)}), '1 2';
is $seq->submeta_text(2, 3), '1 2';
is join (' ', @{$seq->submeta(2, 3, "9 9")}), '9 9';
is $seq->submeta_text(2, 3, "8 8"), '8 8';

is join (' ', @{$seq->named_submeta('trace', 2, 3)}), '5 10';
is $seq->named_submeta_text('trace', 2, 3), '5 10';
is join (' ', @{$seq->named_submeta('trace', 2, 3, "9 9")}), '9 9';
is $seq->named_submeta_text('trace', 2, 3, "8 8"), '8 8';


ok $seq = Bio::Seq::Quality->new(
    -seq => "ATGGGGGTGGTGGTACCCTATGGGGGTGGTGGTACCCT",
    -qual => "10 59 12 75 63 76 84 36 42 10 35 97 81 50 81 53 93 13 38 10 59 12 75 63 76 84 36 42 10 35 97 81 50 81 53 93 13 38",
    -trace_indices => "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38"
);

my $rev;
ok $rev = $seq->revcom; 
is $rev->seq, 'AGGGTACCACCACCCCCATAGGGTACCACCACCCCCAT';
is $rev->qual_text, "38 13 93 53 81 50 81 97 35 10 42 36 84 76 63 75 12 59 10 38 13 93 53 81 50 81 97 35 10 42 36 84 76 63 75 12 59 10";


# selecting ranges based on quality

# test seq with three high quality regions (13, 12 and 3), one very short (3)
ok $seq = Bio::Seq::Quality->new(
    -seq => "ATGGGGGTGGTGGTACCCTATGGGGGTGGTGGTACCCT",
    -qual => "0 5 10 20 30 40 40 50 50 50 50 50 40 10 10 10 5 5 20 20 30 40 50 44 44 50 50 50 50 50 5 5 40 40 40 40 50 50"
    );


is $seq->threshold, undef;
is $seq->threshold(10), 10;
is $seq->threshold(13), 13;

is $seq->count_clear_ranges, 3;

my $newseq = $seq->get_clear_range;
is $newseq->length, 12;


my @ranges = $seq->get_all_clean_ranges;
is scalar @ranges, 3;
my $min_length = 10;
@ranges = $seq->get_all_clean_ranges($min_length);
is scalar @ranges, 2;

my $seqio = Bio::SeqIO->new(
    -file   => test_input_file('test_clear_range.fastq'),
    -format => 'fastq'
);

while ( my $seq = $seqio->next_seq() ) {
    my $newqualobj;
    lives_ok { $newqualobj = $seq->get_clear_range(0) };
    if ($newqualobj) {
        is($newqualobj->id, $seq->id, 'Bug 2845');
    } else {
        ok(0, "No object returned via get_clear_range()");
    }
}




#############################################
#
# try testing some 'meta morphic relations'
#

## belief; As the threshold is increased, the number of clear ranges
## (ncr) should not decrease.

## belief; As the thrshold is increased, the length of the clear
## ranges (lcr) should not decrease.

## belief; As the threshold is incrazed, the clear range length (clr)
## should not increase. Sorry for the terribe var names.

## belief; The number of clear ranges should vary between zero and
## half the sequence length.

## belief; The length of the clear ranges should vary between zero and
## the sequence length.

## belief; The length of the clear range should vary between zero and
## the sequence length.

## belief; The lenght of the clear range should not be larger than the
## length of hte clear ranges.


my @bases = qw (A T C G a t c g);
my @qualities = 0..65;


## See beliefs above:
my $ncr_thresh_sanity = 0;
my $lcr_thresh_sanity = 0;
my $clr_thresh_sanity = 0;

my $ncr_range_sanity = 0;
my $lcr_range_sanity = 0;
my $clr_range_sanity = 0;

my $final_loss_of_sanity = 0;



## Go time:

for (1..100){
    $seq = join("", map {$bases[rand(@bases)]} 1..1000  );
    $qual = join(" ", map {$qualities[rand(@qualities)]} 1..1000  );
    
    $seq = Bio::Seq::Quality->
	new(
	    -seq => $seq,
	    -qual => $qual,
	    );
    
    $seq->threshold(10);
    my $a_ncr = $seq->count_clear_ranges;
    my $a_lcr = $seq->clear_ranges_length;
    my $a_clr = scalar(@{$seq->get_clear_range->qual});
    
    $ncr_range_sanity ++ if $a_ncr >= 0 && $a_ncr <=  500;
    $lcr_range_sanity ++ if $a_lcr >= 0 && $a_lcr <= 1000;
    $clr_range_sanity ++ if $a_clr >= 0 && $a_clr <= 1000;
    $final_loss_of_sanity ++ if $a_lcr >= $a_clr;

    $seq->threshold(20);
    my $b_ncr = $seq->count_clear_ranges;
    my $b_lcr = $seq->clear_ranges_length;
    my $b_clr = scalar(@{$seq->get_clear_range->qual});

    $ncr_range_sanity ++ if $b_ncr >= 0 && $b_ncr <=  500;
    $lcr_range_sanity ++ if $b_lcr >= 0 && $b_lcr <= 1000;
    $clr_range_sanity ++ if $b_clr >= 0 && $b_clr <= 1000;
    $final_loss_of_sanity ++ if $b_lcr >= $b_clr;


    $seq->threshold(30);
    my $c_ncr = $seq->count_clear_ranges;
    my $c_lcr = $seq->clear_ranges_length;
    my $c_clr = scalar(@{$seq->get_clear_range->qual});

    $ncr_range_sanity ++ if $c_ncr >= 0 && $c_ncr <=  500;
    $lcr_range_sanity ++ if $c_lcr >= 0 && $c_lcr <= 1000;
    $clr_range_sanity ++ if $c_clr >= 0 && $c_clr <= 1000;
    $final_loss_of_sanity ++ if $c_lcr >= $c_clr;
    
    
    
    $ncr_thresh_sanity ++ if 
	$a_ncr <= $b_ncr && 
	$b_ncr <= $c_ncr;
    
    $lcr_thresh_sanity ++ if 
	$a_ncr <= $b_ncr && 
	$b_ncr <= $c_ncr;
    
    $clr_thresh_sanity ++ if
	$a_clr >= $b_clr && 
	$b_clr >= $c_clr;
       
}

is $ncr_thresh_sanity, 100;
is $lcr_thresh_sanity, 100;
is $clr_thresh_sanity, 100;

is $ncr_range_sanity, 300;
is $lcr_range_sanity, 300;
is $clr_range_sanity, 300;

is $final_loss_of_sanity, 300;




## Test the mask sequence function ...

## Ideally we'd at least test each function with each permutation of constructors.

my $x = Bio::Seq::Quality->
    new( -seq => "aaaattttccccgggg",
	 -qual =>"1 1 1 1 2 2 2 2 1 1 1 1 3 3 3 3");

$x->threshold(1); 

is $x->mask_below_threshold, "aaaattttccccgggg";

$x->threshold(2); 

is $x->mask_below_threshold, "XXXXttttXXXXgggg";

$x->threshold(3); 

is $x->mask_below_threshold, "XXXXXXXXXXXXgggg";

$x->threshold(4); 

is $x->mask_below_threshold, "XXXXXXXXXXXXXXXX";

