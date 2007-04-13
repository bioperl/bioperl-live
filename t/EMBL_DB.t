# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use vars qw($NUMTESTS $DEBUG);

my $error;

BEGIN { 
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;
	$DEBUG = $ENV{'BIOPERLDEBUG'} || 0;
    $NUMTESTS = 16;
    eval { require IO::String;
			require HTTP::Request::Common;   };
    if( $@ ) {
	    plan skip_all => "IO::String not installed. This means the Bio::DB::* modules are not usable. Skipping tests";
	}
    elsif (!$DEBUG) {
		plan skip_all => 'Must set BIOPERLDEBUG=1 for network tests';
	}
	else {
		plan tests => $NUMTESTS;
	}
	use_ok('Bio::DB::EMBL');
}

my $verbose = 0;

## End of black magic.

my ($db,$seq,$seqio);
# get a single seq

$seq = $seqio = undef;

SKIP: { 
    ok defined($db = new Bio::DB::EMBL(-verbose=>$verbose)); 
    ok(defined($seq = $db->get_Seq_by_acc('J00522')));
    is( $seq->length, 408); 
    ok defined ($db->request_format('fasta'));
	
    eval {ok(defined($seq = $db->get_Seq_by_acc('J02231')))};
	skip('could not connect to embl',2) if $@;
    is( $seq->id, 'embl|J02231|J02231');
    is( $seq->length, 200); 
    ok( defined($db = new Bio::DB::EMBL(-verbose=>$verbose, 
					-retrievaltype => 'tempfile')));
    eval {ok(defined($seqio = $db->get_Stream_by_id(['BUM'])))};
	skip('could not connect to embl',2) if $@;
    undef $db; # testing to see if we can remove gb
    ok( defined($seq = $seqio->next_seq()));
    is( $seq->length, 200);
};

$seq = $seqio = undef;

SKIP: {
    $db = new Bio::DB::EMBL(-verbose => $verbose,
			    -retrievaltype => 'tempfile',
			    -format => 'fasta'
			    ); 
    eval{ok( defined($seqio = $db->get_Stream_by_acc(['J00522 AF303112 J02231'])))};
	skip('could not connect to embl',3) if $@;
    my %seqs;
    # don't assume anything about the order of the sequences
    while ( my $s = $seqio->next_seq ) {
		my ($type,$x,$name) = split(/\|/,$s->display_id);
		$seqs{$x} = $s->length;
    }
    is($seqs{'J00522'},408);
    is($seqs{'AF303112'},1611);
    is($seqs{'J02231'},200);
};
