# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use Test;
use strict;
BEGIN { plan tests => 29 }

use Bio::Variation::Allele;
use Bio::Variation::RNAChange;
use Bio::Variation::DNAMutation;
use Bio::Variation::AAChange;
 
my $obj = Bio::Variation::RNAChange -> new;

ok $obj;

$obj->start(4);           
ok $obj->start, 4;

$obj->end(4); 
ok $obj->end, 4;

$obj->length(1);

ok $obj->length, 1;

$obj->strand('1');  
ok $obj->strand, '1';

ok ($obj->primary_tag, 'Variation' );

$obj->source_tag('source');
ok ($obj->source_tag, 'source' );

$obj->frame(2);   
ok ($obj->frame, 2 );

$obj->score(2);   
ok ($obj->score, 2 );

#test gff string
#$obj->dna_mut('dna_mut'); 
#if ($obj->dna_mut eq 'dna_mut' ) {
#    print "ok 11\n";  
#} else {
#    print "not ok 11\n";
#} 
ok(1);

my $a1 = Bio::Variation::Allele->new(-seq => 'g');
$obj->allele_ori($a1);

ok( $obj->allele_ori->seq, 'g' );

my $a2 = Bio::Variation::Allele->new('-seq' => 'a');
$obj->allele_mut($a2);

ok ($obj->allele_mut->seq, 'a' );

$obj->upStreamSeq('gaagattcagccaagctcaaggatg'); 
ok ($obj->upStreamSeq, 'gaagattcagccaagctcaaggatg' );

$obj->cds_end(1000); 
ok ($obj->cds_end, 1000 );

$obj->dnStreamSeq('aagtgcagttagggctgggaagggt'); 
ok ($obj->dnStreamSeq, 'aagtgcagttagggctgggaagggt' );

$obj->codon_pos(1); 
ok ($obj->codon_pos, 1 );

my $obj3 = Bio::Variation::AAChange -> new;
$obj3->start(2);
$obj->AAChange($obj3);

ok ($obj->label, 'missense' , "label is". $obj->label);


$obj->status('proven'); 
ok ($obj->status, 'proven' );

$obj->proof('experimental'); 
ok ($obj->proof, 'experimental' );
ok ($obj->restriction_changes, '-BccI' );

$obj->region('coding'); 
ok ($obj->region, 'coding' );
$obj->numbering('coding'); 
ok ($obj->numbering, 'coding' );

ok ($obj->codon_ori, 'gaa', "Codon_ori is |". $obj->codon_ori. "|");

ok ($obj->codon_mut, 'aaa' , "Codon_mut is |". $obj->codon_mut. "|");


$obj->codon_pos(1); 
ok ($obj->codon_pos, 1 );
ok( $obj->codon_table, 1 );

$obj->codon_table(3);
ok ( $obj->codon_table, 3 );

$obj->mut_number(2);
ok ( $obj->mut_number, 2 );

$obj->verbose(2);
ok ( $obj->verbose, 2 );
