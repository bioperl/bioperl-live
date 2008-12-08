# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use Data::Dumper;
BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 185);
	
    use_ok('Bio::Tools::Signalp::ExtendedSignalp');
}

###############################################
### TESTS ON SUMMARY OUTPUT FORMAT (NN+HMM) ###
###############################################

my $res = {
	   '1' => {
	           'id'     => 'BC1G_00003.1',
                   'pred'   => 'Signal peptide',
		   'nnpred' => 'signal-peptide',
                   'end'    => '22',
                   'prob'   => '0.999',
                   'anchor' => '0.000',
		  },
        '2' => {
		   'id'     => 'BC1G_00008.1',
                   'pred'   => 'Non-secretory protein',
		   'nnpred' => 'signal-peptide',
                   'end'    => '83',
                   'prob'   => '0.222',
                   'anchor' => '0.067',
		  },
	  };

# Test on filtered results
my $facts   = [qw(maxS D)];
my $in      = test_input_file("signalp.summary");
my $signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
							-file => $in,
							-factors => $facts,
						       );

ok($signalp);
my $i = 1;

while(my $feat = $signalp->next_feature()){
    #print Dumper($feat);
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('anchorProb'))[0], $res->{$i}->{anchor});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}

# Tests without filters.
# It should by default only parses results with Ymax and meanS to mimic default behavior
# from Bio::Tools::Signalp
$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'pred'   => 'Signal peptide',
		'end'    => '22',
		'prob'   => '0.999',
		'anchor' => '0.000',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'pred'   => 'Non-secretory protein',
		'end'    => '83',
		'prob'   => '0.222',
		'anchor' => '0.067',
	       },
       };

#No filters required
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in
						    );
ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('anchorProb'))[0], $res->{$i}->{anchor});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    $i++;
}

#############################################
### TESTS ON SHORT OUTPUT FORMAT (NN+HMM) ###
#############################################

$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'pred'   => 'Signal peptide',
		'nnpred' => 'signal-peptide',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'nnpred' => 'signal-peptide',
		'pred'   => 'Non-secretory protein',
		'end'    => '83',
	       },
	'3' => {
		'id'     => 'BC1G_00009.1',
		'nnpred' => 'signal-peptide',
		'pred'   => 'Non-secretory protein',
		'end'    => '28',
		},
	'4' => {
		'id'     => 'BC1G_00010.1',
		'nnpred' => 'signal-peptide',
		'pred'   => 'Non-secretory protein',
		'end'    => '15',
		},

       };

# Test on filtered results
$facts   = [qw(maxC)];
$in      = test_input_file("signalp.short");
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in,
						     -factors => $facts,
						    );

ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}

# Tests without filters.
# It should by default only parses results with Ymax and meanS to mimic default behavior
# from GPI::Bio::Tools::Signalp
$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'pred'   => 'Signal peptide',
        'prob' => 0.999,
		'nnpred' => 'signal-peptide',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'pred'   => 'Non-secretory protein',
        'prob' => 0.222,
		'nnpred' => 'signal-peptide',
		'end'    => '83',
	       },
       };

#No filters required
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in
						    );
ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('anchorProb'))[0], $res->{$i}->{anchor});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}


###########################################
### TESTS ON SUMMARY OUTPUT FORMAT (NN) ###
###########################################

$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'nnpred' => 'signal-peptide',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'nnpred' => 'signal-peptide',
		'end'    => '83',
	       },
       };

# Test on filtered results BROKEN
$facts   = [qw(maxC)];
$in      = test_input_file("signalp.nn.summary");
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in,
						     -factors => $facts,
						    );

ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){

    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}
# Tests without filters.
# It should by default only parses results with Ymax and meanS to mimic default behavior
# from GPI::Bio::Tools::Signalp
$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'nnpred' => 'signal-peptide',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'nnpred' => 'signal-peptide',
		'end'    => '83',
	       },
	'3' => {
		'id'     => 'BC1G_00009.1',
		'nnpred' => 'signal-peptide',
		'end'    => '28',
		},
	'4' => {
		'id'     => 'BC1G_00010.1',
		'nnpred' => 'signal-peptide',
		'end'    => '15',
		},
       };

#No filters required BROKEN
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in
						    );
ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){

    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}


############################################
### TESTS ON SUMMARY OUTPUT FORMAT (HMM) ###
############################################

$res = {
	'1' => {
		'id'     => 'BC1G_00002.1',
		'prob'   => '0.000',
		'anchor' => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00003.1',
		'prob'   => '0.999',
		'anchor' => '0.000',
		'cleav'  => '0.973',
		'pred'   => 'Signal peptide',
		'end'    => '22',
	       },
	'3' => {
		'id'     => 'BC1G_00004.1',
		'prob'   => '0.003',
		'anchor' => '0.000',
		'cleav'  => '0.001',
		'pred'   => 'Non-secretory protein',
		'end'    => '19',
	       },
	'4' => {
		'id'     => 'BC1G_00005.1',
		'prob'   => '0.008',
		'anchor' => '0.000',
		'cleav'  => '0.007',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'5' => {
		'id'     => 'BC1G_00006.1',
		'prob'   => '0.000',
		'anchor' => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '23',
		},
	'6' => {
		'id'     => 'BC1G_00007.1',
		'prob'   => '0.240',
		'anchor' => '0.000',
		'cleav'  => '0.228',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'7' => {
		'id'     => 'BC1G_00008.1',
		'prob'   => '0.222',
		'anchor' => '0.067',
		'cleav'  => '0.061',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'8' => {
		'id'     => 'BC1G_00009.1',
		'prob'   => '0.000',
		'anchor' => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '20',
		},
       };

# It is impossible to filter with hmm output...
$in      = test_input_file("signalp.hmm.summary");
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in,
						    );

ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('anchorProb'))[0], $res->{$i}->{anchor});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    is(($feat->get_tag_values('cleavageSiteProb'))[0], $res->{$i}->{cleav});
    $i++;
}

#########################################
### TESTS ON SHORT OUTPUT FORMAT (NN) ###
#########################################

$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'nnpred' => 'signal-peptide',
		'yprob'  => '0.866',
		'dprob'  => '0.902',
		'cprob'  => '0.934',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'yprob'  => '0.383',
		'dprob'  => '0.436',
		'cprob'  => '0.576',
		'nnpred' => 'signal-peptide',
		'end'    => '83',
	       },
       };

# Test on filtered results
$facts   = [qw(maxY)];
$in      = test_input_file("signalp.nn.short");
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in,
						     -factors => $facts,
						    );

ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('maxCprob'))[0], $res->{$i}->{cprob});
    is(($feat->get_tag_values('Dprob'))[0], $res->{$i}->{dprob});
    is(($feat->get_tag_values('maxYprob'))[0], $res->{$i}->{yprob});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}
# Tests without filters.
# It should by default only parses results with Ymax and meanS to mimic default behavior
# from GPI::Bio::Tools::Signalp
$res = {
	'1' => {
		'id'     => 'BC1G_00003.1',
		'nnpred' => 'signal-peptide',
		'yprob'  => '0.866',
		'dprob'  => '0.902',
		'cprob'  => '0.934',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00008.1',
		'yprob'  => '0.383',
		'dprob'  => '0.436',
		'cprob'  => '0.576',
		'nnpred' => 'signal-peptide',
		'end'    => '83',
	       },
       };

#No filters required
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in
						    );
ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('maxCprob'))[0], $res->{$i}->{cprob});
    is(($feat->get_tag_values('Dprob'))[0], $res->{$i}->{dprob});
    is(($feat->get_tag_values('maxYprob'))[0], $res->{$i}->{yprob});
    is(($feat->get_tag_values('nnPrediction'))[0], $res->{$i}->{nnpred});
    $i++;
}

##########################################
### TESTS ON SHORT OUTPUT FORMAT (HMM) ###
##########################################

$res = {
	'1' => {
		'id'     => 'BC1G_00002.1',
		'prob'   => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
	       },
	'2' => {
		'id'     => 'BC1G_00003.1',
		'prob'   => '0.999',
		'cleav'  => '0.973',
		'pred'   => 'Signal peptide',
		'end'    => '22',
	       },
	'3' => {
		'id'     => 'BC1G_00004.1',
		'prob'   => '0.003',
		'cleav'  => '0.001',
		'pred'   => 'Non-secretory protein',
		'end'    => '19',
	       },
	'4' => {
		'id'     => 'BC1G_00005.1',
		'prob'   => '0.008',
		'cleav'  => '0.007',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'5' => {
		'id'     => 'BC1G_00006.1',
		'prob'   => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '23',
		},
	'6' => {
		'id'     => 'BC1G_00007.1',
		'prob'   => '0.240',
		'cleav'  => '0.228',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'7' => {
		'id'     => 'BC1G_00008.1',
		'prob'   => '0.222',
		'cleav'  => '0.061',
		'pred'   => 'Non-secretory protein',
		'end'    => '22',
		},
	'8' => {
		'id'     => 'BC1G_00009.1',
		'prob'   => '0.000',
		'cleav'  => '0.000',
		'pred'   => 'Non-secretory protein',
		'end'    => '20',
		},
       };

# No filters available with hmm on short output
$in      = test_input_file("signalp.hmm.short");
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in,
						    );

ok($signalp);
$i = 1;

while(my $feat = $signalp->next_feature()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('cleavageSiteProb'))[0], $res->{$i}->{cleav});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    $i++;
}


exit 0;
