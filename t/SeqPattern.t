# -*-Perl-*-
## Bioperl Test Harness Script for Modules
##

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
    plan tests => 6;
}

use Bio::Tools::SeqPattern;

my ( $pattern,$pattern_obj,$pattern_obj2, $pattern_obj3);

$pattern     = '(CCCCT)N{1,200}(agyyg)N{1,80}(ag)';
$pattern_obj = new Bio::Tools::SeqPattern(-SEQ =>$pattern, -TYPE =>'dna');
ok defined($pattern_obj) && ref($pattern_obj) && $pattern_obj->isa('Bio::Tools::SeqPattern');

$pattern_obj2  = $pattern_obj->revcom();
ok $pattern_obj2->str, '(CT)N(CRRCT){1,80}N(AGGGG){1,200}';

$pattern_obj3 = $pattern_obj->revcom(1);
ok $pattern_obj3->str, '(CT).{1,80}(C[GA][GA]CT).(AGGGG){1,200}';

$pattern     = '(CCCCT)N{1,200}(agyyg)N{1,80}(bb)'; # test protein object expand
$pattern_obj = new Bio::Tools::SeqPattern(-SEQ =>$pattern, -TYPE =>'protein');
ok defined($pattern_obj) && ref($pattern_obj) && $pattern_obj->isa('Bio::Tools::SeqPattern');

ok $pattern_obj2->expand, '(CT).(C[AG][AG]CT){1,80}.(AGGGG){1,200}';

# amino patterns

$pattern = 'ABZH';
$pattern_obj2 = new Bio::Tools::SeqPattern(-SEQ =>$pattern, 
					   -TYPE =>'amino');
ok $pattern_obj2->expand, 'A[EQ][DN]H';
