# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use ExtUtils::testlib;
use strict;
require 'dumpvar.pl';

BEGIN {
	# to handle systems with no installed Test module
	# we include the t dir (where a copy of Test.pm is located)
	# as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 18;
}

        # redirect STDERR to STDOUT
open (STDERR, ">&STDOUT");

print("Checking if the Bio::Seq::SeqWithQuality module could be used...\n");
        # test 1
use Bio::Seq::SeqWithQuality;
ok(1);

use Bio::PrimarySeq;
use Bio::Seq::PrimaryQual;

# create some random sequence object with no id
my $seqobj_broken = Bio::PrimarySeq->new( -seq => "ATCGATCGA",
                            );
	# dumpValue($seqobj_broken);

my $seqobj = Bio::PrimarySeq->new( -seq => "ATCGATCGA",
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                            );
ok(!$@);




# create some random quality object with the same number of qualities and the same identifiers
my $string_quals = "10 20 30 40 50 40 30 20 10";
my $qualobj;
eval {
$qualobj = Bio::Seq::PrimaryQual->new( -qual => $string_quals,
                            -id  => 'QualityFragment-12',
                            -accession_number => 'X78121',
                            );
};
ok(!$@);

# check to see what happens when you construct the SeqWithQuality object
my $swq1 = Bio::Seq::SeqWithQuality->new( -seq	=>	$seqobj,
					-qual		=>	$qualobj);
ok(!$@);

print("Testing various weird constructors...\n");
print("\ta) No ids, Sequence object, no quality...\n");
	# w for weird
my $wswq1;
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq  =>	$seqobj,
						-qual	=>	"");
};
ok(!$@);


print("\tb) No ids, no sequence, quality object...\n");
	# note that you must provide a alphabet for this one.
$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
					-qual => $qualobj,
					-alphabet => 'dna'
);
print("\tc) Absolutely nothing. (HAHAHAHA)...\n");
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
						-qual => "",
						-alphabet => 'dna'
	);
};
ok(!$@);
print("\td) Absolutely nothing but an ID\n");
eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq => "",
						-qual => "",
						-alphabet => 'dna',
						-id => 'an object with no sequence and no quality but with an id'
	);
};
ok(!$@);

print("\td) No sequence, No quality, No ID...\n");

eval {
	$wswq1 = Bio::Seq::SeqWithQuality->new( -seq  =>	"",
							-qual	=>	"");
};
	# this should fail without a alphabet
ok($@);
	# dumpValue($wswq1);





print("Testing various methods and behaviors...\n");

print("1. Testing the seq() method...\n");
	print("\t1a) get\n");
	my $original_seq = $swq1->seq();
	ok ($original_seq eq "ATCGATCGA");
	print("\t1b) set\n");
	ok ($swq1->seq("AAAAAAAAAAAA"));
	print("\t1c) get (again, to make sure the set was done.)\n");
	ok($swq1->seq() eq "AAAAAAAAAAAA");
	print("\tSetting the sequence back to the original value...\n");
	$swq1->seq($original_seq);

print("2. Testing the qual() method...\n");
	print("\t2a) get\n");
	my @qual = @{$swq1->qual()};
	my $str_qual = join(' ',@qual);
	ok ($str_qual eq "10 20 30 40 50 40 30 20 10");
	print("\t2b) set\n");
	ok ($swq1->qual("10 10 10 10 10"));
	print("\t2c) get (again, to make sure the set was done.)\n");
	my @qual2 = @{$swq1->qual()};
	my $str_qual2 = join(' ',@qual2);
	ok($str_qual2 eq "10 10 10 10 10");
	print("\tSetting the quality back to the original value...\n");
	$swq1->qual($str_qual);

print("3. Testing the length() method...\n");
	print("\t3a) When lengths are equal...\n");
	ok($swq1->length() == 9);	
	print("\t3b) When lengths are differenct\n");
	$swq1->qual("10 10 10 10 10");
	ok($swq1->length() eq "DIFFERENT");

print("4. Testing the qual_obj() method...\n");
	print("\t4a) Testing qual_obj()...\n");
		my $retr_qual_obj = $swq1->qual_obj();
		ok (ref($retr_qual_obj) eq "Bio::Seq::PrimaryQual");
	print("\t4b) Testing qual_obj(\$ref)...\n");
		$swq1->qual_obj($qualobj);

print("5. Testing the seq_obj() method...\n");
	print("\t5a) Testing seq_qual_obj()...\n");
		my $retr_seq_obj = $swq1->seq_obj();
		ok (ref($retr_seq_obj) eq "Bio::PrimarySeq");
	print("\t5b) Testing seq_obj(\$ref)...\n");
		$swq1->seq_obj($seqobj);



