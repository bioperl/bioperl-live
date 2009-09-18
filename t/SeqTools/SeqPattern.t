# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 28);
	
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

SKIP: {
    test_skip(-tests => 19, -requires_module => 'List::MoreUtils');
    # Test reverse complement
    my $rev_pattern = $pattern_obj2->backtranslate;
    isa_ok $rev_pattern, 'Bio::Tools::SeqPattern';
    is $rev_pattern->str, 'GCNRAYSARCAY';

    # Test exceptions.
    throws_ok { $pattern_obj2->revcom  } qr/revcom for .+ sequence types/;
    throws_ok { $rev_pattern->backtranslate } qr/backtranslate for .+ sequence types/;

    # Test reverse translation more thoroughly
    
    my @data;
    while (<DATA>) {
       chomp;
       push @data, [ split ];
    }
    
    foreach my $line (@data) {
        my $pattern_obj = Bio::Tools::SeqPattern->new(
            -SEQ  => $line->[0],
            -TYPE => 'amino',
        );
    
        isa_ok $pattern_obj, 'Bio::Tools::SeqPattern';
    
        my $backtranslate = $pattern_obj->backtranslate;
        isa_ok $backtranslate, 'Bio::Tools::SeqPattern';
        is $backtranslate->str, $line->[1],
    }
}

__DATA__
MAEELKAVAP ATGGCNGARGARYTNAARGCNGTNGCNCCN
LKGHB[WhYq]Q YTNAARGGNCAYRAYYRNCAR
(LK){2,3}[^GHB][WHYQ]Q (YTNAAR){2,3}HBNYRNCAR
LK[^GHB][WHYQ]Q YTNAARHBNYRNCAR
(LK){2,3}[^GHB][WHYQ]QX.X (YTNAAR){2,3}HBNYRNCARNNNNNNNNN
