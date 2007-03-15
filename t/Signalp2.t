# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Signalp2.t,v 1.1 cjfields 2007/03/14 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use Data::Dumper;
my $DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;
	plan tests => 32;
    use_ok('Bio::Tools::Signalp::ExtendedSignalp');
    use_ok('Bio::Root::IO');
}

my $res = {
	   '1' => {
	           'id' => 'BC1G_00003.1',
                   'pred' => 'Signal peptide',
                   'end'  => '21',
                   'prob' => '1.000',
                   'anchor' => '0.000',
		  },
           '2' => {
		   'id' => 'BC1G_00008.1',
                   'pred' => 'Signal anchor',
                   'end'  => '82',
                   'prob' => '0.194',
                   'anchor' => '0.737',
		  },
	  };

# Test on filtered results

my $facts   = [qw(maxS D)];
my $in      = Bio::Root::IO->catfile("t","data","exsignalp.out");
my $signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
							-file => $in,
							-factors => $facts,
						       );

ok($signalp);
my $i = 1;

while(my $feat = $signalp->next_result()){
    is($feat->seq_id(), $res->{$i}->{id});
    is($feat->end(), $res->{$i}->{end});
    is(($feat->get_tag_values('peptideProb'))[0], $res->{$i}->{prob});
    is(($feat->get_tag_values('anchorProb'))[0], $res->{$i}->{anchor});
    is(($feat->get_tag_values('signalpPrediction'))[0], $res->{$i}->{pred});
    $i++;
}

# Tests without filters.
$res = {
	'1'  => { 'id' => 'BC1G_00001.1' , 'maxC_score' => '0.074' },
	'2'  => { 'id' => 'BC1G_00002.1' , 'maxC_score' => '0.079' },
	'3'  => { 'id' => 'BC1G_00003.1' , 'maxC_score' => '0.934' },
	'4'  => { 'id' => 'BC1G_00004.1' , 'maxC_score' => '0.078' },
	'5'  => { 'id' => 'BC1G_00005.1' , 'maxC_score' => '0.085' },
	'6'  => { 'id' => 'BC1G_00006.1' , 'maxC_score' => '0.244' },
	'7'  => { 'id' => 'BC1G_00007.1' , 'maxC_score' => '0.257' },
	'8'  => { 'id' => 'BC1G_00008.1' , 'maxC_score' => '0.576' },
	'9'  => { 'id' => 'BC1G_00009.1' , 'maxC_score' => '0.425' },
	'10' => { 'id' => 'BC1G_00010.1' , 'maxC_score' => '0.054' },
       };

#No filters required
$signalp = Bio::Tools::Signalp::ExtendedSignalp->new(
						     -file => $in
						    );
ok($signalp);
$i = 1;

while(my $feat = $signalp->next_result()){
    is($feat->seq_id(), $res->{$i}->{id});
    is(($feat->get_tag_values('maxC_score'))[0], $res->{$i}->{maxC_score});
    $i++;
}

exit 0;
