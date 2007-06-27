# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib 't/lib';
    use BioperlTest;
    
    test_begin(-tests => 9);
	
	use_ok('Bio::Tools::SeqPattern');
}

my ( $pattern,$pattern_obj,$pattern_obj2, $pattern_obj3);

$pattern     = '(CCCCT)N{1,200}(agyyg)N{1,80}(ag)';
ok $pattern_obj = Bio::Tools::SeqPattern->new(-SEQ =>$pattern, -TYPE =>'dna');
isa_ok $pattern_obj, 'Bio::Tools::SeqPattern';

$pattern_obj2  = $pattern_obj->revcom();
is $pattern_obj2->str, '(CT)N(CRRCT){1,80}N(AGGGG){1,200}';

$pattern_obj3 = $pattern_obj->revcom(1);
is $pattern_obj3->str, '(CT).{1,80}(C[GA][GA]CT).(AGGGG){1,200}';

$pattern     = '(CCCCT)N{1,200}(agyyg)N{1,80}(bb)'; # test protein object expand
ok $pattern_obj = Bio::Tools::SeqPattern->new(-SEQ =>$pattern, -TYPE =>'protein');
isa_ok $pattern_obj, 'Bio::Tools::SeqPattern';

is $pattern_obj2->expand, '(CT).(C[AG][AG]CT){1,80}.(AGGGG){1,200}';

# amino patterns

$pattern = 'ABZH';
$pattern_obj2 = Bio::Tools::SeqPattern->new(-SEQ =>$pattern, 
					   -TYPE =>'amino');
is $pattern_obj2->expand, 'A[EQ][DN]H';
