# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 53);
	
	use_ok('Bio::Seq::Quality');
}

my $DEBUG = test_debug();

# create some random sequence object with no id
my $seqobj_broken = Bio::Seq::Quality->new( -seq => "ATCGATCGA",
                                          );

my $seqobj;
lives_ok {
	$seqobj = Bio::Seq::Quality->new( -seq => "ATCGATCGA",
                                     -id  => 'QualityFragment-12',
                                     -accession_number => 'X78121',
                                   );
};


# create some random quality object with the same number of qualities and the same identifiers
my $string_quals = "10 20 30 40 50 40 30 20 10";
my $qualobj;
lives_ok {
	$qualobj = Bio::Seq::Quality->new( -qual => $string_quals,
                                      -id  => 'QualityFragment-12',
                                      -accession_number => 'X78121',
                                        );
};

# check to see what happens when you construct the Quality object
ok my $swq1 = Bio::Seq::Quality->new( -seq => "ATCGATCGA",
                                      -id  => 'QualityFragment-12',
                                      -accession_number => 'X78121',
                                      -qual	=>	$string_quals);


print("Testing various weird constructors...\n") if $DEBUG;
print("\ta) No ids, Sequence object, no quality...\n") if $DEBUG;
# w for weird
my $wswq1;
lives_ok {
	$wswq1 = Bio::Seq::Quality->new( -seq  =>  "ATCGATCGA",
                                         -qual	=>	"");
};
print $@ if $DEBUG;


print("\tb) No ids, no sequence, quality object...\n") if $DEBUG;
	# note that you must provide a alphabet for this one.
$wswq1 = Bio::Seq::Quality->new( -seq => "",
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
    $wswq1 = Bio::Seq::Quality->new( -seq => "",
                                            -qual => "",
                                            -alphabet => 'dna',
                                            -id => 'an object with no sequence and no quality but with an id'
	);
};

print("\td) No sequence, No quality, No ID...\n") if $DEBUG;
warnings_like { ok $wswq1 = Bio::Seq::Quality->new( -seq  =>	"",
                                    -qual =>	"",
                                    -verbose => 0
) } qr/Got a sequence with no letters in it cannot guess alphabet/;

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
          # ok ('1 2 3' eq join(' ',@{$swq1->subqual(1,3)}));
     print("\t6d) Testing the subqual at the end (border condition)\n") if $DEBUG;
          # ok ('7 8 9' eq join(' ',@{$swq1->subqual(7,9)}));
     print("\t6d) Testing the subqual in the middle\n") if $DEBUG;
          # ok ('4 5 6' eq join(' ',@{$swq1->subqual(4,6)}));

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
is $swq1->length, $swq2->length;


#
# end of test inherited from seqwithquality.t 
#
#################################################################
#
# testing new functionality
#

my $qual = '0 1 2 3 4 5 6 7 8 9 11 12';
my $trace = '0 5 10 15 20 25 30 35 40 45 50 55';

ok my $seq = Bio::Seq::Quality->new
    ( -qual => $qual,
      -trace_indices => $trace,
      -seq =>  'atcgatcgatcg',
      -id  => 'human_id',
      -accession_number => 'S000012',
      -verbose => $DEBUG >= 0 ? $DEBUG : 0
);

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
