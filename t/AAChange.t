# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'


use Test;
use strict;

BEGIN { plan tests => 24 }
use Bio::Variation::Allele;
use Bio::Variation::AAChange;
use Bio::Variation::RNAChange;


## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 


my $obj = Bio::Variation::AAChange -> new;
ok(1);
ok defined $obj;
ok ref($obj), 'Bio::Variation::AAChange';

$obj->start(3);           
ok $obj->start, 3;

$obj->end(3); 
ok $obj->end, 3;

$obj->length(3);

ok $obj->length, 3;

$obj->strand('1');  
ok $obj->strand, '1';

ok $obj->primary_tag, 'Variation';

$obj->source_tag('source');
ok $obj->source_tag, 'source';

$obj->frame(2);   
ok $obj->frame,2;

$obj->score(2);   
ok $obj->score, 2;

$obj->isMutation(1); 
ok $obj->isMutation;

my $a1 = Bio::Variation::Allele->new(-seq => 'V');
$obj->allele_ori($a1);

ok $obj->allele_ori->seq, 'V';

my $a2 = Bio::Variation::Allele->new('-seq' => 'A');
$obj->add_Allele($a2);

ok $obj->allele_mut->seq, 'A';

$obj->upStreamSeq('upStreamSeq'); 
ok $obj->upStreamSeq, 'upStreamSeq';

$obj->dnStreamSeq('dnStreamSeq'); 
ok $obj->dnStreamSeq, 'dnStreamSeq' ;

ok $obj->label, 'substitution';

$obj->status('proven'); 
ok $obj->status, 'proven';

$obj->proof('experimental'); 
ok $obj->proof, 'experimental';

$obj->region('region'); 
ok $obj->region, 'region';

$obj->region_value('region_value'); 
ok $obj->region_value, 'region_value';

$obj->numbering('coding'); 
ok $obj->numbering, 'coding';

my $obj2 = Bio::Variation::RNAChange -> new(-start => 7, 
					  -end => 7,
					  -cds_end => 100,
					  -codon_pos => 1,
					  -upStreamSeq => 'acgcgcgcgc',
					  -dnStreamSeq => 'acgcgcgcgc'
					  );
$obj2->label('missense');
$obj->RNAChange($obj2);

ok $obj->trivname, 'V3A', "Trivial name is !". $obj->trivname. "!\n";

$obj->mut_number(2);
ok $obj->mut_number, 2;
