# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Domcut.t,v 1.1 2003/07/23 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG $ERROR $METAERROR);
BEGIN {
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    eval { require Test::More; };
    if( $@ ) {
		use lib 't/lib';
    }
    use Test::More;

    $NUMTESTS = 27;
	eval {
		 require IO::String; 
		 require LWP::UserAgent;
	};
	if( $@ ) {
		plan skip_all => "IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests";
	} elsif (!$DEBUG) {
		plan skip_all => 'Must set BIOPERLDEBUG=1 for network tests';
	} else {
	    plan tests => $NUMTESTS;
	}
	use_ok('Bio::Seq::Meta::Array');
	use_ok('Bio::PrimarySeq');
	use_ok('Bio::Tools::Analysis::Protein::Domcut');
}

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);

######## test using PrimarySeq object ##############
my $seq = Bio::PrimarySeq->new(-seq        => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQPPPPPPPPPPPPPDQRS',
							   -display_id => 'test2');

ok $tool = Bio::Tools::Analysis::Protein::Domcut->new( -seq=>$seq);
ok $tool->run ();
if ($tool->status eq 'TERMINATED_BY_ERROR') {
	skip('Problem with DomCut run, check status', 21);
}
ok my $raw    = $tool->result('');
ok my $parsed = $tool->result('parsed');
is ($parsed->[23]{'score'}, '-0.209');
my @res       = $tool->result('Bio::SeqFeatureI');
if (scalar @res > 0) {
	ok 1;
} else {
	skip('No network access - could not connect to Domcut server', 18);
}
ok my $meta = $tool->result('meta');

if (!$METAERROR) { #if Bio::Seq::Meta::Array available
is($meta->named_submeta_text('Domcut', 1,2), "0.068 0.053");
is ($meta->seq, "MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQPPPPPPPPPPPPPDQRS");
}
	
	
########## test using Bio::Seq object ##############
ok my $tool2 = Bio::WebAgent->new(-verbose =>$verbose);

ok my $seq2  = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
			 -display_id => 'test2');

ok $tool2 = Bio::Tools::Analysis::Protein::Domcut->new( -seq=>$seq2->primary_seq);
ok $tool2->run ();

@res = $tool2->result('Bio::SeqFeatureI');

if (scalar @res > 0) {
	ok 1;
} else {
	skip('No network access - could not connect to Domcut server', 10);
}

ok my $parsed2 = $tool2->result('parsed');
is ($parsed2->[23]{'score'}, '-0.209');

ok my $meta2 = $tool2->result('meta');

is($meta2->named_submeta_text('Domcut', 1,2), "0.068 0.053");
is ($meta2->seq, "MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS");

ok my $seq4 = new Bio::Seq;
ok $seq2->primary_seq($meta2);
ok $seq2->add_SeqFeature(@res);
ok $seq2->primary_seq->named_submeta_text('Domcut', 1,2);

