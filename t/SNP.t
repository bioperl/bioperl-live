# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

use strict;

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) { 
	use lib 't\lib';
    }
    use Test::More;
    plan tests => 13;
	use_ok('Bio::Variation::SNP');
}

my($a);

#
# SNP
#

ok $a = new Bio::Variation::SNP;
is $a->id('123'), 123;
eval { $a->di('123'); };
ok 1 if $@;
is $a->validated('by-cluster'), 'by-cluster';
my @alleles = ('A', 'T');
is $a->validated(\@alleles), \@alleles;
is $a->desc('abc'), 'abc'; # Bio::Variation::Allele method
is $a->chromosome('X'), 'X'; # Bio::Variation::Allele method
ok my $s = $a->add_subsnp;
ok $s->is_subsnp;
is $s->handle('HGBASE'), 'HGBASE';
ok $a->add_subsnp;
is $a->each_subsnp, 2;


