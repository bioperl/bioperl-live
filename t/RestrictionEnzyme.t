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
    plan test => 14 }

use Bio::Tools::RestrictionEnzyme;
use Bio::Seq;

ok(1);

# pre declare 6cutters to prevent warnings

my (@sixcutters,$seq,$dna,$l1,$l2,@fragments,$re);

$dna = 'CCTCCGGGGACTGCCGTGCCGGGCGGGAATTCGCCATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGA';

# Build sequence and restriction enzyme objects.
ok defined( $seq = new Bio::Seq(-ID  =>'test_seq',
			     -SEQ =>$dna)); 

ok $seq->id();
ok $re  = new Bio::Tools::RestrictionEnzyme(-NAME=>'EcoRI');
ok $re->seq->id;
ok $re->site;
ok $re->palindromic;
ok scalar(@fragments = $re->cut_seq($seq)), 2, "can't cut sequence";
ok scalar @fragments, 2;

ok $l1 = length $fragments[0], 27;
ok $l2 = length $fragments[1], 43;
ok $l1 + $l2, $seq->length, "sum of frag lengths != length of original seq\n";

ok @sixcutters = $re->available_list(6), 304, "can't get list of 6-cutters";
ok $re->is_available('HindIII');




