# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#
# modeled after the t/Allele.t test script

use ExtUtils::testlib;
use strict;
use Dumpvalue qw(dumpValue);


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

my $dumper = new Dumpvalue();

print("Checking to see if Bio::Tools::Primer3 is available.\n");
use Bio::Tools::Primer3;
ok(1);

print("Checking to see if Bio::SeqFeature::Primer is available.\n");
use Bio::SeqFeature::Primer;
ok(1);

print("Checking to see if a primer outfile can be read...\n");
     # ph = "primer handle"
my $ph = Bio::Tools::Primer3->new('-file' => Bio::Root::IO->catfile("t","data",
          "primer3_outfile.txt"));
ok(ref($ph) eq "Bio::Tools::Primer3");
print("Getting an object from the primer file...\n");

     # now get a primer!
my $thingy = $ph->next_primer();
print("That thingy isa $thingy\n");
ok (ref($thingy) eq "Bio::Seq::PrimedSeq");
# print("Here is the primer object that I got back:\n");
# $dumper->dumpValue($thingy);
# 
# sub feature_things {
#         my @features = $thingy->all_SeqFeatures();
#         print("These are the seqfeatures: @features\n");
#         print("Dumping out those features names...\n");
#         foreach (@features) {
#              print("Name: ".$_->seqname()."\n");
#              print("\tSequence ".$_->entire_seq()->seq()."\n");
#         }
# }
# 



sub the_old_module {
}
