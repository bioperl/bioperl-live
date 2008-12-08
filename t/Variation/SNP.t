# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13);
	
	use_ok('Bio::Variation::SNP');
}

my($a);

#
# SNP
#

ok $a = Bio::Variation::SNP->new();
is $a->id('123'), 123;
eval { $a->di('123'); };
ok $@;
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
