# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------

use strict;
use Test;

BEGIN { plan tests => 39}

use Bio::Variation::SeqDiff;
use Bio::Variation::DNAMutation;
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

ok $obj->dna_ori('gctgctgatcgatcgtagctagctag'); 
ok $obj->dna_ori, 'gctgctgatcgatcgtagctagctag';

ok $obj->dna_mut('gctgctgatcggtcgtagctagctag'); 
ok $obj->dna_mut, 'gctgctgatcggtcgtagctagctag';

ok $obj->rna_ori('gctgctgatcgatcgtagctagctag'); 
ok $obj->rna_ori, 'gctgctgatcgatcgtagctagctag';

$obj->rna_mut('gctgctgatcgatcgtagctagctag'); 
ok $obj->rna_mut, 'gctgctgatcgatcgtagctagctag';

ok $obj->aa_ori('MHYTRD'); 
ok $obj->aa_ori, 'MHYTRD';

ok $obj->aa_mut('MHGTRD'); 
ok $obj->aa_mut, 'MHGTRD';

ok $m = Bio::Variation::DNAMutation->new;
ok $obj->add_Variant($m);

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
