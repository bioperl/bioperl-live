# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 8);
	
    use_ok('Bio::SeqFeature::Primer');
}

my $DEBUG = test_debug();

print("Checking to see if a BSFP object can be created:\n") if $DEBUG;
# yes sure, but first scope a few variables
my $seqsequence = "gcatcgatctagctagcta";
my $primersequence = "aaaaaacgatcgatcgtagctagct";

my $seqname = "chads_nifty_sequence";
my $primername = "chads_nifty_primer";
# ok, and what about variables governing where the feature is located?
# check the primer3docs, luke...
# TARGET=513,26
# PRIMER_FIRST_BASE_INDEX=1
# PRIMER_LEFT=484,20


print("Checking to see if the BSFP object can be constructed with a bio::seq object\n") if $DEBUG;
my $seq = Bio::Seq->new( -seq => $seqsequence, -id =>$seqname);
my $bsfp_seq = Bio::SeqFeature::Primer->new( -sequence => $seq,
                                             -TARGET => '5,3' );
isa_ok $bsfp_seq, "Bio::SeqFeature::Primer";

print("Checking to see if the BSFP object can be constructed with scalars\n") if $DEBUG;

my $bsfp_scalar = Bio::SeqFeature::Primer->new( -sequence => $primersequence,
                                        -id => $primername,
                                             -TARGET => '5,3' );
isa_ok $bsfp_scalar, "Bio::SeqFeature::Primer";

print("Checking to see that seq() returns a Bio::Seq object and that the object is the right one.\n") if $DEBUG;
isa_ok $bsfp_scalar->seq(), "Bio::Seq";
print("First for the scalar-ily created one.\n") if $DEBUG;
print("id ok?\n") if $DEBUG;
is $bsfp_scalar->seq()->id(), $primername;
print("sequence ok?\n") if $DEBUG;
is $bsfp_scalar->seq()->seq(), $primersequence;
print("Now for the seq-ily created one\n") if $DEBUG;
print("id ok?\n") if $DEBUG;
is $bsfp_seq->seq()->display_id(), $seqname;
print("sequence ok?\n") if $DEBUG;
is $bsfp_seq->seq()->seq(), $seqsequence;
