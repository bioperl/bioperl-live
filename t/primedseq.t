# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use strict;
#use Dumpvalue qw(dumpValue);
use vars qw($DEBUG);


BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 4;
}

$DEBUG = $ENV{'BIOPERLDEBUG'};
#my $dumper = new Dumpvalue();

use Bio::Seq;
use Bio::SeqFeature::Primer;

print("Checking to see if Bio::Seq::PrimedSeq is available.\n") if $DEBUG;
use Bio::Seq::PrimedSeq;
ok(1);
print("Trying to create a PrimedSeq...\n") if $DEBUG;
print("First create a seq object.\n") if $DEBUG;
my $seqobj = new Bio::Seq( -seq => "aaaaaaaaaaaaactgatcgatcgatcg",
                         -display_name => "chads_kewl_sequence");
ok (ref($seqobj) eq "Bio::Seq");
print("Now create a left primer...\n") if $DEBUG;
my $left_primer = new Bio::SeqFeature::Primer( -sequence => "tttttttttttagctgtgca",
                                             -id => "left_primer"
                                        );
ok (ref($left_primer) eq "Bio::SeqFeature::Primer");
print("Now create a right primer...\n") if $DEBUG;
my $right_primer = new Bio::SeqFeature::Primer( -sequence => "gggggggggggggcacgtcgat",
                                             -id => "right_primer");
ok (ref($right_primer) eq "Bio::SeqFeature::Primer");
print("Now create the primedseq object...\n") if $DEBUG;
my $ps = new Bio::Seq::PrimedSeq( -target_sequence => $seqobj,
                                   -left_primer => $left_primer,
                                   -right_primer => $right_primer);
print("isa $ps\n") if $DEBUG;
# ok (ref($ps) eq "Bio::Seq::PrimedSeq");




