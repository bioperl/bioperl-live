# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

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
    plan tests => 42;
}

use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
use Bio::Variation::Allele;

ok 1;
my ($obj, $mm, $aa, $dna, $m);

ok $obj = Bio::Variation::SeqDiff -> new;

ok $obj->id('id');
ok $obj->id, 'id';

ok $obj->sysname('sysname');
ok $obj->sysname, 'sysname';

$obj->trivname('trivname'); 
ok $obj->trivname eq 'trivname';

ok $obj->chromosome('chr');
ok $obj->chromosome, 'chr';

ok $obj->description('desc');
ok $obj->description, 'desc';

ok $obj->numbering('numbering');
ok $obj->numbering, 'numbering';

ok $obj->offset(100);
ok $obj->offset, 100;
#                  12345678901234567890
ok $obj->dna_ori ('gctgctgatcgatcgtagctagctag');
ok $obj->dna_ori, 'gctgctgatcgatcgtagctagctag';

# generate mutated DNA seq from the mutation
ok $m = Bio::Variation::DNAMutation->new(-isMutation => 1, -start=>14, -end=>14);
ok $a = Bio::Variation::Allele->new(-seq=>'c');
$b = Bio::Variation::Allele->new(-seq=>'g');
ok $m->allele_ori($a);
ok $m->allele_mut($b);
ok $obj->add_Variant($m);
my $m2 = Bio::Variation::DNAMutation->new(-isMutation => 1, -start=>19, -end=>19);
my $a2 = Bio::Variation::Allele->new(-seq=>'c');
my $b2 = Bio::Variation::Allele->new(-seq=>'g');
$m2->allele_ori($a2);
$m2->allele_mut($b2);
$obj->add_Variant($m2);

#ok $obj->dna_mut('gctgctgatcggtcgtagctagctag');
ok $obj->dna_mut, 'gctgctgatcgatggtaggtagctag';

ok $obj->rna_ori('gctgctgatcgatcgtagctagctag');
ok $obj->rna_ori, 'gctgctgatcgatcgtagctagctag';

$obj->rna_mut('gctgctgatcgatcgtagctagctag'); 
ok $obj->rna_mut, 'gctgctgatcgatcgtagctagctag';

ok $obj->aa_ori('MHYTRD');
ok $obj->aa_ori, 'MHYTRD';

ok $obj->aa_mut('MHGTRD');
ok $obj->aa_mut, 'MHGTRD';

foreach $mm ( $obj->each_Variant ) {
    $mm->primary_tag('a');
    ok $mm->isa('Bio::Variation::VariantI');
}


ok $obj->gene_symbol('fos');
ok $obj->gene_symbol, 'fos';

ok $obj->rna_offset(10);
ok $obj->rna_offset == 10;

ok $obj->rna_id('transcript#3');
ok $obj->rna_id, 'transcript#3';

ok $dna = $obj->seqobj('dna_ori');
ok $dna->isa('Bio::PrimarySeq');

$obj->aa_mut(''); 
$aa = $obj->seqobj('aa_mut');
ok not defined $aa;

eval {
    $dna = $obj->seqobj('dna_ri');
};
ok $@;
