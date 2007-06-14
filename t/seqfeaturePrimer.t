# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use strict;
#use Dumpvalue qw(dumpValue);

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

my $DEBUG = $ENV{'BIOPERLDEBUG'};
#my $dumper = new Dumpvalue();

print("Checking to see if Bio::SeqFeature::Primer is available.\n") if $DEBUG;
use Bio::SeqFeature::Primer;
ok(1);
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
ok(ref($bsfp_seq) eq "Bio::SeqFeature::Primer");

print("Checking to see if the BSFP object can be constructed with scalars\n") if $DEBUG;

my $bsfp_scalar = Bio::SeqFeature::Primer->new( -sequence => $primersequence,
                                        -id => $primername,
                                             -TARGET => '5,3' );
ok(ref($bsfp_scalar) eq "Bio::SeqFeature::Primer");

print("Checking to see that seq() returns a Bio::Seq object and that the object is the right one.\n") if $DEBUG;
ok(ref($bsfp_scalar->seq()) eq "Bio::Seq");
print("First for the scalar-ily created one.\n") if $DEBUG;
print("id ok?\n") if $DEBUG;
ok($bsfp_scalar->seq()->id() eq $primername);
print("sequence ok?\n") if $DEBUG;
ok($bsfp_scalar->seq()->seq() eq $primersequence);
print("Now for the seq-ily created one\n") if $DEBUG;
print("id ok?\n") if $DEBUG;
ok($bsfp_seq->seq()->display_id() eq $seqname);
print("sequence ok?\n") if $DEBUG;
ok($bsfp_seq->seq()->seq() eq $seqsequence);

print("Here is the structure of the BSFP_scalar object:\n") if $DEBUG;
# $dumper->dumpValue($bsfp_scalar);


