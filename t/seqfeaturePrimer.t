# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use ExtUtils::testlib;
use strict;
use Dumpvalue qw(dumpValue);
use Bio::SeqIO;


BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 8;
}

my $dumper = new Dumpvalue();

print("Checking to see if Bio::SeqFeature::Primer is available.\n");
use Bio::SeqFeature::Primer;
ok(1);
print("Checking to see if a BSFP object can be created:\n");
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


print("Checking to see if the BSFP object can be constructed with a bio::seq object\n");
my $seq = new Bio::Seq( -seq => $seqsequence, -display_id =>$seqname);
my $bsfp_seq = new Bio::SeqFeature::Primer( -sequence => $seq,
                                             -TARGET => '5,3' );
ok(ref($bsfp_seq) eq "Bio::SeqFeature::Primer");

print("Checking to see if the BSFP object can be constructed with scalars\n");
     # 

my $bsfp_scalar = new Bio::SeqFeature::Primer( -sequence => $primersequence,
                                        -id => $primername,
                                             -TARGET => '5,3' );
ok(ref($bsfp_scalar) eq "Bio::SeqFeature::Primer");

print("Checking to see that seq() returns a Bio::Seq object and that the object is the right one.\n");
ok(ref($bsfp_scalar->seq()) eq "Bio::Seq");
print("First for the scalar-ily created one.\n");
ok($bsfp_scalar->seq()->display_id() eq $primername);
ok($bsfp_scalar->seq()->seq() eq $primersequence);
print("Now for the seq-ily created one\n");
ok($bsfp_seq->seq()->display_id() eq $seqname);
ok($bsfp_seq->seq()->seq() eq $seqsequence);

print("Here is the structure of the BSFP_scalar object:\n");
$dumper->dumpValue($bsfp_scalar);


