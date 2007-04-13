# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: MitoProt.t,v 1.1 2003/07/26 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG $ERROR);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
	eval { require Test::More; };
	$ERROR = 0;
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;

	$NUMTESTS = 10;
	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if( $@ ) {
		plan skip_all => "IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests";
	} elsif (!$DEBUG) {
        plan skip_all => 'Must set BIOPERLDEBUG=1 for network tests';
    } else {
		plan tests => $NUMTESTS
	};
	use_ok 'Bio::Tools::Analysis::Protein::Mitoprot';
	use_ok 'Bio::PrimarySeq';
	use_ok 'Bio::WebAgent';
}

use Data::Dumper;

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seq = Bio::PrimarySeq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDSFFGSDFDGDS'.
                               'DFGSDFGSDGDFGSDFGDSFGDGFSDRSRQDQRS',
                               -display_id => 'test2');

ok $tool = Bio::Tools::Analysis::Protein::Mitoprot->new( -seq=>$seq);
ok $tool->run ();
exit if $tool->status eq 'TERMINATED_BY_ERROR';
ok my $raw = $tool->result('');
ok my $parsed = $tool->result('parsed');
is ($parsed->{'charge'}, -13);
ok my @res = $tool->result('Bio::SeqFeatureI');

