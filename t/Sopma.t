# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Sopma.t,v 1.1 2003/07/23 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	$NUMTESTS = 16;
	
	eval { require Test::More; };
	if( $@ ) {
		use lib 't/lib';
	}
	use Test::More;

	plan tests => $NUMTESTS;

	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if ($@) {
		plan skip_all => "IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests";
	}
	
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::Analysis::Protein::Sopma');
}

my $verbose = $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

my $seq = Bio::PrimarySeq->new(
  -seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
  -display_id => 'test2');
ok $tool = Bio::Tools::Analysis::Protein::Sopma->new( -seq=>$seq,
                                                      -window_width => 15);

SKIP: {
	skip "Skipping tests which require remote servers - set env variable BIOPERLDEBUG to test", 12 unless $DEBUG;
	
	ok $tool->run();
	skip "Tool was terminated by some error: problem connecting to server?", 11 if $tool->status eq 'TERMINATED_BY_ERROR';
	
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is ($parsed->[0]{'helix'}, '104');
	ok my @res = $tool->result('Bio::SeqFeatureI');
	ok my $meta = $tool->result('meta', "ww15");

	ok $tool->window_width(21);
	ok $tool->clear();
	ok $tool->run;
	ok my $meta2 = $tool->result('meta', "ww21");
	
	SKIP: {
		eval {
			require Bio::Seq::Meta::Array;
		};
		skip "Bio::Seq::Meta::Array not installed - will skip tests using meta sequences", 2 if $@;
		is $meta->named_submeta_text('Sopma_helix|ww15',1,2), '104 195';
		is $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS';
	}
}
