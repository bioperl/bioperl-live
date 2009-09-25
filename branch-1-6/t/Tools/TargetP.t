# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 124);
	
	use_ok('Bio::Tools::TargetP');
}

my $targetp = Bio::Tools::TargetP->new(
				       -file => test_input_file('targetp.out')
				      );

ok($targetp);

my $items = {
	     '1' => {
		     'id' => 'BC1G_00001.1',
		     'len' => '173',
		     'mTP' => '0.393',
		     'SP'  => '0.024',
		     'other' => '0.683',
		     'loc'   => '_',
		     'RC'    => '4',
		     'tplen' => undef,
		    },
	     '2' => {
		     'id' => 'BC1G_00002.1',
		     'len' => '120',
		     'mTP' => '0.619',
		     'SP'  => '0.040',
		     'other' => '0.458',
		     'loc'   => 'M',
		     'RC'    => '5',
		     'tplen' => '97',
		    },
	     '3' => {
		     'id' => 'BC1G_00003.1',
		     'len' => '323',
		     'mTP' => '0.094',
		     'SP'  => '0.895',
		     'other' => '0.027',
		     'loc'   => 'S',
		     'RC'    => '1',
		     'tplen' => '21',
		    },
	     '4' => {
		     'id' => 'BC1G_00004.1',
		     'len' => '361',
		     'mTP' => '0.402',
		     'SP'  => '0.072',
		     'other' => '0.479',
		     'loc'   => '_',
		     'RC'    => '5',
		     'tplen' => undef,
		    },
	     '5' => {
		     'id' => 'BC1G_00005.1',
		     'len' => '244',
		     'mTP' => '0.526',
		     'SP'  => '0.035',
		     'other' => '0.548',
		     'loc'   => '_',
		     'RC'    => '5',
		     'tplen' => undef,
		    },
	     '6' => {
		     'id' => 'BC1G_00006.1',
		     'len' => '35',
		     'mTP' => '0.234',
		     'SP'  => '0.036',
		     'other' => '0.819',
		     'loc'   => '_',
		     'RC'    => '3',
		     'tplen' => undef,
		    },
	     '7' => {
		     'id' => 'BC1G_00007.1',
		     'len' => '73',
		     'mTP' => '0.292',
		     'SP'  => '0.127',
		     'other' => '0.431',
		     'loc'   => '_',
		     'RC'    => '5',
		     'tplen' => undef,
		    },
	     '8' => {
		     'id' => 'BC1G_00008.1',
		     'len' => '349',
		     'mTP' => '0.088',
		     'SP'  => '0.958',
		     'other' => '0.078',
		     'loc'   => 'S',
		     'RC'    => '1',
		     'tplen' => '82',
		    },
	     '9' => {
		     'id' => 'BC1G_00009.1',
		     'len' => '514',
		     'mTP' => '0.183',
		     'SP'  => '0.102',
		     'other' => '0.735',
		     'loc'   => '_',
		     'RC'    => '3',
		     'tplen' => undef,
		    },
	     '10' => {
		     'id' => 'BC1G_00010.1',
		     'len' => '440',
		     'mTP' => '0.114',
		     'SP'  => '0.088',
		     'other' => '0.865',
		     'loc'   => '_',
		     'RC'    => '2',
		     'tplen' => undef,
		    },
	     '11' => {
		     'id' => 'BC1G_04501.1',
		     'len' => '215',
		     'mTP' => '0.185',
		     'SP'  => '0.038',
		     'other' => '0.843',
		     'loc'   => '_',
		     'RC'    => '2',
		     'tplen' => undef,
		    },
	     '12' => {
		     'id' => 'BC1G_04502.1',
		     'len' => '395',
		     'mTP' => '0.118',
		     'SP'  => '0.164',
		     'other' => '0.825',
		     'loc'   => '_',
		     'RC'    => '2',
		     'tplen' => undef,
		    },
	     '13' => {
		     'id' => 'BC1G_04503.1',
		     'len' => '199',
		     'mTP' => '0.515',
		     'SP'  => '0.062',
		     'other' => '0.436',
		     'loc'   => 'M',
		     'RC'    => '5',
		     'tplen' => '20',
		    },
	     '14' => {
		     'id' => 'BC1G_04504.1',
		     'len' => '220',
		     'mTP' => '0.440',
		     'SP'  => '0.030',
		     'other' => '0.707',
		     'loc'   => '_',
		     'RC'    => '4',
		     'tplen' => undef,
		    },
	     '15' => {
		     'id' => 'BC1G_04505.1',
		     'len' => '67',
		     'mTP' => '0.382',
		     'SP'  => '0.049',
		     'other' => '0.610',
		     'loc'   => '_',
		     'RC'    => '4',
		     'tplen' => undef,
		    },
	    };
	
my $i = 1;
$targetp->_parse_results();

is($targetp->network(), 'NON-PLANT');
is($targetp->analysis_method_version(), "v1.1");

while(my $feat = $targetp->next_prediction()){

    is($feat->seq_id(), $items->{$i}->{id}, "good SeqID");
    is($feat->length(), $items->{$i}->{len}, "good Seqlength");
    is(($feat->get_tag_values('mitochondrionCutOff'))[0], $items->{$i}->{mTP}, "correct Mitochondrion cutoff");
    is(($feat->get_tag_values('signalPeptideCutOff'))[0], $items->{$i}->{SP}, "correct signalpPeptide cutoff");
    is(($feat->get_tag_values('otherCutOff'))[0], $items->{$i}->{other}, "correct other cutoff");
    is(($feat->get_tag_values('location'))[0], $targetp->_toString_location($items->{$i}->{loc}), "correct location");
    is(($feat->get_tag_values('reliabilityClass'))[0], $items->{$i}->{RC}, "correct Reliability class score");
    
	if ($feat->has_tag('signalPeptideLength')) {
		is(($feat->get_tag_values('signalPeptideLength'))[0], $items->{$i}->{tplen}, "correct peptide signal length")
	} else {
		is($feat->has_tag('signalPeptideLength'), '', "No peptide signal length reported")
	}
    $i++;
}
