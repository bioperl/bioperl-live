# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

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
    plan tests => 3;
}


print("Checking if the Bio::SeqIO::Qual module could be used, even though it shouldn't be directly use'd...\n");
        # test 1
use Bio::SeqIO::qual;
ok(1);

print("Checking to see if PrimaryQual.pm can be used...\n");
use Bio::Seq::PrimaryQual;
ok(1);

print("Checking to see if PrimaryQual objects can be created from a file...\n");
my $in_qual  = Bio::SeqIO->new(-file => "<t/data/qualfile.qual" , '-format' => 'qual');
ok(1);

my @quals;
print("I saw these in qualfile.qual:\n");
while ( my $qual = $in_qual->next_qual() ) {
		# ::dumpValue($qual);
	print($qual->id()."\n");
	@quals = @{$qual->qual()};
	print("(".scalar(@quals).") quality values.\n");
}
