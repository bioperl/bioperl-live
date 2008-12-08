# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 31);
	
	use_ok('Bio::Variation::Allele');
	use_ok('Bio::Variation::RNAChange');
	use_ok('Bio::Variation::AAChange');
}

ok my $obj = Bio::Variation::RNAChange->new();

$obj->start(4);           
is $obj->start, 4;

$obj->end(4); 
is $obj->end, 4;

$obj->length(1);

is $obj->length, 1;

$obj->strand('1');  
is $obj->strand, '1';

is ($obj->primary_tag, 'Variation' );

$obj->source_tag('source');
is ($obj->source_tag, 'source' );

$obj->frame(2);   
is ($obj->frame, 2 );

$obj->score(2);   
is ($obj->score, 2 );

#test gff string
#$obj->dna_mut('dna_mut'); 
#if ($obj->dna_mut eq 'dna_mut' ) {
#    print "ok 11\n";  
#} else {
#    print "not ok 11\n";
#}

my $a1 = Bio::Variation::Allele->new(-seq => 'g');
$obj->allele_ori($a1);

is( $obj->allele_ori->seq, 'g' );

my $a2 = Bio::Variation::Allele->new('-seq' => 'a');
$obj->allele_mut($a2);

is($obj->allele_mut->seq, 'a' );

$obj->upStreamSeq('gaagattcagccaagctcaaggatg'); 
is ($obj->upStreamSeq, 'gaagattcagccaagctcaaggatg' );

$obj->cds_end(1000); 
is ($obj->cds_end, 1000 );

$obj->dnStreamSeq('aagtgcagttagggctgggaagggt'); 
is ($obj->dnStreamSeq, 'aagtgcagttagggctgggaagggt' );

$obj->codon_pos(1); 
is ($obj->codon_pos, 1 );

my $obj3 = Bio::Variation::AAChange->new();
$obj3->start(2);
$obj->AAChange($obj3);
$obj3->allele_ori($a1);
$obj3->allele_mut($a2);

is ($obj->label, 'missense' , "label is". $obj->label);


$obj->status('proven'); 
is ($obj->status, 'proven' );

$obj->proof('experimental'); 
is ($obj->proof, 'experimental' );
is ($obj->restriction_changes, '-BccI' );

$obj->region('coding'); 
is ($obj->region, 'coding' );
$obj->numbering('coding'); 
is ($obj->numbering, 'coding' );

is ($obj->codon_ori, 'gaa', "Codon_ori is |". $obj->codon_ori. "|");

is($obj->codon_mut, 'aaa' , "Codon_mut is |". $obj->codon_mut. "|");


$obj->codon_pos(1); 
is ($obj->codon_pos, 1 );
is( $obj->codon_table, 1 );

$obj->codon_table(3);
is ( $obj->codon_table, 3 );

$obj->mut_number(2);
is ( $obj->mut_number, 2 );

$obj->verbose(2);
is ( $obj->verbose, 2 );
