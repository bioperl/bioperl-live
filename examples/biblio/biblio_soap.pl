#!/usr/bin/perl
#
# This was actually a part of the test suite - but because it starts
# an external process it was safer not to use it as a test (the process
# could be left running if an error occurs).
#
# It is an example of a TCP-based SOAP exchange.
#

use strict;
eval { require SOAP::Lite;
};
if( $@ ){
    die("must have SOAP::Lite installed to run this script");
}

use vars qw($NUMTESTS);

my $error;

BEGIN { 
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    $error = 0;
    if( $@ ) {
	use lib 't';
    }
    use Test;
    plan tests => 10;
}

my $testnum;
my $verbose = 0;

use Bio::Biblio;

# --- launch a testing SOAP server
my ($pid, $port, $max_port);
$port = 4444;
$max_port = $port + 100;
if ($pid = fork) {
    # parent here

    sleep 1;
    my $biblio = new Bio::Biblio (-location => "tcp://localhost:$port",
				  -namespace => 'soap_server');
    
    ok ($biblio->get_count, '43');
    ok ($biblio->get_by_id ('X'), 'X');
    ok ($biblio->find ('a,b','c,d')->get_collection_id, 'a,b,c,d');
    ok ($biblio->find (['x', 'y'], ['u', 'v'])->get_collection_id, 'x,y,u,v');

    ok ( eval { join (',', @{ $biblio->find ('AAA')->get_all_ids }) }, 'AAA'); print STDERR $@ if $@;

    ok ( eval { join (',', @{ $biblio->find ('XXX')->get_all }) }, 'XXX'); print STDERR $@ if $@;

    ok ( eval { $biblio->find (46)->has_next }, 1); print STDERR $@ if $@;

    ok ( eval { $biblio->find ('BBB')->get_next }, 'BBB'); print STDERR $@ if $@;

    ok ( eval { join (',', @{ $biblio->find ('CCC')->get_more (3) }) }, 'CCC,CCC,CCC'); print STDERR $@ if $@;

    ok ( eval { $biblio->find (46)->exists }, 0); print STDERR $@ if $@;


    # clean-up the running server
    kill 9, $pid if defined $pid;
    print "    SOAP server $pid killed\n";

} elsif (defined $pid) {
    # child here - a testing SOAP server

    package soap_server;
    use strict;
    use SOAP::Transport::TCP;
    my $daemon;
    while ($port < $max_port) {
	eval {
	    $daemon = SOAP::Transport::TCP::Server
		-> new (LocalAddr => 'localhost', LocalPort => $port, Listen => 5, Reuse => 1)
		    -> dispatch_to('soap_server');
	};
	last unless $@;
	$port++;
    }
    print "    Contact to SOAP server at ", join(':', $daemon->sockhost, $daemon->sockport), " (server PID: $$)\n";
    $daemon->handle;

    sub getBibRefCount { shift;  return 43; }
    sub getById { shift; return shift; }
    sub find {
	my ($self, $keywords, $attrs) = @_;
	return join (',', (@{ $keywords }, @{ $attrs })) if $attrs;
	return join (',', @{ $keywords });
    }
    sub getAllIDs { shift; return [ shift ] }
    sub getAllBibRefs { shift; return [ shift ] }
    sub hasNext { return SOAP::Data->type (boolean => 'true'); }
    sub getNext { shift; return [ '1', shift]; }
    sub getMore {
	my ($self, $id, $how_many) = @_;
	my @result = ('1');
	push (@result, $id) for (1..$how_many);
	return \@result;
    }
    sub exists { return SOAP::Data->type (boolean => '0'); }
    sub destroy {}

    package main;

} else {
        # fork failed
        print STDERR "Testing SOAP services FAILED: $!.\n";
    }
