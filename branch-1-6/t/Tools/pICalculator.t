# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 38);
	
    use_ok('Bio::Seq');
    use_ok('Bio::Tools::pICalculator');
}

my @results = (12.999052267583,12.99700393539,12.9905348815881,12.9701609055248,12.9065486239062,12.7131376670492,12.1681721433832,10.8960154975975,8.82939162036317,6.81329734996812,5.58311842185452,4.87361913724596,4.11053952923425,3.00644711484741,1.91237900622079,1.19755236429121,0.669596284738213,0.0571988207175853,-0.584285455699191,-1.14218959353989,-1.79865831607402,-2.74360327055112,-3.87361697725167,-4.91494976791445,-6.01005299841696,-7.43711791135299,-8.77859455006782,-9.53905973773058,-9.84470802408586);

my $protein = "MVLLLILSVLLLKEDVRGSAQSSERRVVAHMPGDIIIGALFSVHHQPTVDKVHERKCGAVREQYGI";
ok my $seq = Bio::Seq->new(-seq => $protein);
is $seq->seq, $protein;
ok my $pep = $seq->seq;
ok my $calc = Bio::Tools::pICalculator->new(-places => 2);
ok $calc->seq($seq);
ok my $iep = $calc->iep;
for ( my $x = 0 ; $x <= 14 ; $x += .5 ) {
   float_is($calc->charge_at_pH($x), $results[(2 * $x)]);
}
is ($calc->iep,8.54);
