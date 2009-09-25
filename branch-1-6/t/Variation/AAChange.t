# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 27);
	
	use_ok('Bio::Variation::Allele');
	use_ok('Bio::Variation::AAChange');
	use_ok('Bio::Variation::RNAChange');
}


ok my $obj = Bio::Variation::AAChange->new();
isa_ok $obj, 'Bio::Variation::AAChange';

$obj->start(3);           
is $obj->start, 3;

$obj->end(3); 
is $obj->end, 3;

$obj->length(3);

is $obj->length, 3;

$obj->strand('1');  
is $obj->strand, '1';

is $obj->primary_tag, 'Variation';

$obj->source_tag('source');
is $obj->source_tag, 'source';

$obj->frame(2);   
is $obj->frame,2;

$obj->score(2);   
is $obj->score, 2;

$obj->isMutation(1); 
ok $obj->isMutation;

my $a1 = Bio::Variation::Allele->new(-seq => 'V');
$obj->allele_ori($a1);

is $obj->allele_ori->seq, 'V';

my $a2 = Bio::Variation::Allele->new('-seq' => 'A');
$obj->add_Allele($a2);

is $obj->allele_mut->seq, 'A';

is $obj->similarity_score, 0;

$obj->upStreamSeq('upStreamSeq'); 
is $obj->upStreamSeq, 'upStreamSeq';

$obj->dnStreamSeq('dnStreamSeq'); 
is $obj->dnStreamSeq, 'dnStreamSeq' ;

is $obj->label, 'substitution, conservative';

$obj->status('proven'); 
is $obj->status, 'proven';

$obj->proof('experimental'); 
is $obj->proof, 'experimental';

$obj->region('region'); 
is $obj->region, 'region';

$obj->region_value('region_value'); 
is $obj->region_value, 'region_value';

$obj->numbering('coding'); 
is $obj->numbering, 'coding';

my $obj2 = Bio::Variation::RNAChange->new(-start => 7, 
					  -end => 7,
					  -cds_end => 100,
					  -codon_pos => 1,
					  -upStreamSeq => 'acgcgcgcgc',
					  -dnStreamSeq => 'acgcgcgcgc'
					  );
$obj2->label('missense');
$obj->RNAChange($obj2);

is $obj->trivname, 'V3A', "Trivial name is [". $obj->trivname. "]";

$obj->mut_number(2);
is $obj->mut_number, 2;
