# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) { 
	use lib 't';
    }
    use Test;
    plan tests => 13 }


use Bio::Variation::SNP;

ok(1);

my($a);

#
# SNP
#

ok $a = new Bio::Variation::SNP;
ok $a->id('123'), 123;
eval { $a->di('123'); };
ok 1 if $@;
ok $a->validated('by-cluster'), 'by-cluster';
my @alleles = ('A', 'T');
ok $a->validated(\@alleles), \@alleles;
ok $a->desc('abc'), 'abc'; # Bio::Variation::Allele method
ok $a->chromosome('X'), 'X'; # Bio::Variation::Allele method
ok my $s = $a->add_subsnp;
ok $s->is_subsnp;
ok $s->handle('HGBASE'), 'HGBASE';
ok $a->add_subsnp;
ok $a->each_subsnp, 2;


