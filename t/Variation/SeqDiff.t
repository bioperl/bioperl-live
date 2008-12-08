# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 44);
	
	use_ok('Bio::Variation::SeqDiff');
	use_ok('Bio::Variation::DNAMutation');
	use_ok('Bio::Variation::Allele');
}

my ($obj, $mm, $aa, $dna, $m);

ok $obj = Bio::Variation::SeqDiff->new();

ok $obj->id('id');
is $obj->id, 'id';

ok $obj->sysname('sysname');
is $obj->sysname, 'sysname';

$obj->trivname('trivname'); 
is $obj->trivname, 'trivname';

ok $obj->chromosome('chr');
is $obj->chromosome, 'chr';

ok $obj->description('desc');
is $obj->description, 'desc';

ok $obj->numbering('numbering');
is $obj->numbering, 'numbering';

ok $obj->offset(100);
is $obj->offset, 100;
#                  12345678901234567890
ok $obj->dna_ori ('gctgctgatcgatcgtagctagctag');
is $obj->dna_ori, 'gctgctgatcgatcgtagctagctag';

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
is $obj->dna_mut, 'gctgctgatcgatggtaggtagctag';

ok $obj->rna_ori('gctgctgatcgatcgtagctagctag');
is $obj->rna_ori, 'gctgctgatcgatcgtagctagctag';

$obj->rna_mut('gctgctgatcgatcgtagctagctag'); 
is $obj->rna_mut, 'gctgctgatcgatcgtagctagctag';

ok $obj->aa_ori('MHYTRD');
is $obj->aa_ori, 'MHYTRD';

ok $obj->aa_mut('MHGTRD');
is $obj->aa_mut, 'MHGTRD';

foreach $mm ( $obj->each_Variant ) {
    $mm->primary_tag('a');
    isa_ok($mm,'Bio::Variation::VariantI');
}


ok $obj->gene_symbol('fos');
is $obj->gene_symbol, 'fos';

ok $obj->rna_offset(10);
is $obj->rna_offset, 10;

ok $obj->rna_id('transcript#3');
is $obj->rna_id, 'transcript#3';

ok $dna = $obj->seqobj('dna_ori');
isa_ok($dna,'Bio::PrimarySeq');

$obj->aa_mut(''); 
$aa = $obj->seqobj('aa_mut');
ok not defined $aa;

eval {
    $dna = $obj->seqobj('dna_ri');
};
ok $@;
