# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id: Sopma.t,v 1.1 2003/07/23 
# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'
use strict;
use vars qw($NUMTESTS $DEBUG $ERROR $METAERROR);
$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $ERROR = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;

    $NUMTESTS = 12;
    plan tests => $NUMTESTS;

    eval {
	require IO::String; 
	require LWP::UserAgent;
	
    }; 
    if( $@ ) {
        warn("IO::String or LWP::UserAgent not installed. This means that the module is not usable. Skipping tests");
	$ERROR = 1;
    }
	#check this is available, set error flag if not.
	eval {
		require Bio::Seq::Meta::Array;
		};
	if ($@) {
		warn ("Bio::Seq::Meta::Array not installed - will skip tests using meta sequences");
		$METAERROR = 1;
		}
}

END {
    foreach ( $Test::ntest..$NUMTESTS) {
	skip('unable to run all of the tests depending on web access',1);
    }
}

exit 0 if $ERROR ==  1;

use Data::Dumper;

use Bio::Seq;
require Bio::Tools::Analysis::Protein::Sopma;

ok 1;

my $verbose = 0;
$verbose = 1 if $DEBUG;

ok my $tool = Bio::WebAgent->new(-verbose =>$verbose);



my $seq = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
						-display_id => 'test2',
						);
ok $tool = Bio::Tools::Analysis::Protein::Sopma->new( -seq=>$seq,
													 -window_width => 15);
ok $tool->run ();
ok my $raw = $tool->result('');
print "$raw\n";
ok my $parsed = $tool->result('parsed');
ok ($parsed->[0]{'helix'}, '104');
my @res = $tool->result('Bio::SeqFeatureI');
if (scalar @res > 0) {
    ok 1;
} else {
    skip('No network access - could not connect to NetPhos server', 1);
}
ok my $meta = $tool->result('all');

if (!$METAERROR) { #if Bio::Seq::Meta::Array available
	#meta sequence contains data...
	ok ($meta->{'primary_seq'}{'seq'}, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS');

	#but not available thru method call...??
	skip ("meta sequences are undefined?", $meta->named_submeta_text('Sopma_helix',1,2), '104 195');
	skip ("meta sequences are undefined?", $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS');
	}
