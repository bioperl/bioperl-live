# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 22);
	
	use_ok('Bio::Seq::SeqWithQuality');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Seq::PrimaryQual');
}

my $DEBUG = test_debug();

my $verbosity = $DEBUG || -1;

# create some random sequence object with no id
my $seqobj_broken = Bio::PrimarySeq->new( -seq => "ATCGATCGA");

ok my $seqobj = Bio::PrimarySeq->new( -seq => "ATCGATCGA",
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                            -verbose => $verbosity);

# create some random quality object with the same number of qualities and the same identifiers
my $string_quals = "10 20 30 40 50 40 30 20 10";
my $indices = "5 10 15 20 25 30 35 40 45";
my $qualobj;
eval {
$qualobj = Bio::Seq::PrimaryQual->new( -qual => $string_quals,
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                            -verbose => $verbosity);
};
ok(!$@);


# check to see what happens when you construct the SeqWithQuality object
my $swq1 = Bio::Seq::SeqWithQuality->new( -seq	=>	$seqobj,
                                         -verbose => $verbosity,
					-qual		=>	$qualobj);
ok(!$@);
no warnings;

print("Testing various weird constructors...\n") if $DEBUG;
print("\ta) No ids, Sequence object, no quality...\n") if $DEBUG;
# w for weird
my $wswq1;
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq  =>	$seqobj,
                                                -verbose => $verbosity,
						-qual	=>	"");
};
ok(!$@);

print("\tb) No ids, no sequence, quality object...\n") if $DEBUG;
	# note that you must provide a alphabet for this one.
$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
                                        -verbose => $verbosity,
					-qual => $qualobj,
					-alphabet => 'dna'
);
print("\tc) Absolutely nothing. (HAHAHAHA)...\n") if $DEBUG;
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
                                                -verbose => $verbosity,
						-qual => "",
						-alphabet => 'dna'
	);
};
ok(!$@);
print("\td) Absolutely nothing but an ID\n") if $DEBUG;
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
                                                -verbose => $verbosity,
						-qual => "",
						-alphabet => 'dna',
						-id => 'an object with no sequence and no quality but with an id'
	);
};
ok(!$@);

print("\td) No sequence, No quality, No ID...\n") if $DEBUG;

eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq  =>	"",
                                                -verbose => $verbosity,
							-qual	=>	"");
};
# this should fail without a alphabet
ok($@);

print("Testing various methods and behaviors...\n") if $DEBUG;

print("1. Testing the seq() method...\n") if $DEBUG;
	print("\t1a) get\n") if $DEBUG;
	my $original_seq = $swq1->seq();
	is ($original_seq, "ATCGATCGA");
	print("\t1b) set\n") if $DEBUG;
	ok ($swq1->seq("AAAAAAAAAAAA"));
	print("\t1c) get (again, to make sure the set was done.)\n") if $DEBUG;
	is ($swq1->seq(), "AAAAAAAAAAAA");
	print("\tSetting the sequence back to the original value...\n") if $DEBUG;
	$swq1->seq($original_seq);

print("2. Testing the qual() method...\n") if $DEBUG;
	print("\t2a) get\n") if $DEBUG;
	my @qual = @{$swq1->qual()};
	my $str_qual = join(' ',@qual);
	is ($str_qual, "10 20 30 40 50 40 30 20 10");
	print("\t2b) set\n") if $DEBUG;
	ok ($swq1->qual("10 10 10 10 10"));
	print("\t2c) get (again, to make sure the set was done.)\n") if $DEBUG;
	my @qual2 = @{$swq1->qual()};
	my $str_qual2 = join(' ',@qual2);
	is($str_qual2, "10 10 10 10 10");
	print("\tSetting the quality back to the original value...\n") if $DEBUG;
	$swq1->qual($str_qual);

print("3. Testing the length() method...\n") if $DEBUG;
	print("\t3a) When lengths are equal...\n") if $DEBUG;
	is($swq1->length(), 9);	
	print("\t3b) When lengths are different\n") if $DEBUG;
	$swq1->qual("10 10 10 10 10");
	is($swq1->length(), "DIFFERENT");


print("4. Testing the qual_obj() method...\n") if $DEBUG;
	print("\t4a) Testing qual_obj()...\n") if $DEBUG;
		my $retr_qual_obj = $swq1->qual_obj();
		isa_ok $retr_qual_obj, "Bio::Seq::PrimaryQual";
	print("\t4b) Testing qual_obj(\$ref)...\n") if $DEBUG;
		$swq1->qual_obj($qualobj);

print("5. Testing the seq_obj() method...\n") if $DEBUG;
	print("\t5a) Testing seq_qual_obj()...\n") if $DEBUG;
		my $retr_seq_obj = $swq1->seq_obj();
		isa_ok $retr_seq_obj, "Bio::PrimarySeq";
	print("\t5b) Testing seq_obj(\$ref)...\n") if $DEBUG;
		$swq1->seq_obj($seqobj);

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
$swq1 = Bio::Seq::SeqWithQuality->new(-seq =>  'G',
                                      -qual => '0',
                                      -verbose => $verbosity,
                                     );
my $swq2 = Bio::Seq::SeqWithQuality->new(-seq =>  'G',
                                         -qual => '65',
                                         -verbose => $verbosity,
                                     );
is $swq1->length, $swq2->length;

$swq1 = Bio::Seq::SeqWithQuality->new(-seq =>  'GC',
                                      -verbose => $verbosity,
                                      -qual => '0 0',
                                     );
$swq2 = Bio::Seq::SeqWithQuality->new(-seq =>  'GT',
                                      -verbose => $verbosity,
                                      -qual => '65 0',
                                     );
is $swq1->length, $swq2->length;
