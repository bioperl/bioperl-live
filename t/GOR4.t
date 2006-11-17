# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: GOR4.t,v 1.1 2003/07/23 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG);

BEGIN {
	$NUMTESTS = 14;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
	
	eval {require Test::More;};
	if ($@) {
		use lib 't/lib';
	}
	use Test::More;
	
	eval {
		require IO::String; 
		require LWP::UserAgent;
	};
	if ($@) {
		plan skip_all => 'IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	
	use_ok("Data::Dumper");
	use_ok("Bio::Seq");
	use_ok("Bio::Tools::Analysis::Protein::GOR4");
}

#	eval {require Bio::Seq::Meta::Array;};
#	"Bio::Seq::Meta::Array not installed - will skip tests using meta sequences"

my $verbose = $DEBUG;

my $seq = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
                        -display_id => 'test2');
ok my $tool = Bio::Tools::Analysis::Protein::GOR4->new(-seq=>$seq->primary_seq);

SKIP: {
	skip "Skipping tests which require remote servers, set BIOPERLDEBUG=1 to test", 10 unless $DEBUG;
    ok $tool->run();
	skip "Skipping tests since we got terminated by a server error", 9 if $tool->status eq 'TERMINATED_BY_ERROR';
    ok my $raw = $tool->result('');
    ok my $parsed = $tool->result('parsed');
    is $parsed->[0]{'coil'}, '999';
    my @res = sort {$a->start <=> $b->start} $tool->result('Bio::SeqFeatureI');
    if (scalar @res > 0) {
		ok 1;
    }
	else {
		skip 'No results - could not connect to GOR4 server?', 6;
    }
	is $res[0]->start, 1;
	is $res[0]->end, 43;
    ok my $meta = $tool->result('meta');
    
    eval {require Bio::Seq::Meta::Array;};
	skip "Bio::Seq::Meta::Array not installed - will skip tests using meta sequences", 2 if $@;
	is $meta->named_submeta_text('GOR4_coil',1,2), '999 999';
	is $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS';
}
