# This is -*-Perl-*- code
## Bioperl Test Harness Script for Modules
##
# $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Biblio.t'

use strict;
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

    plan tests => 25;
}

my $testnum;
my $verbose = 0;

## End of black magic.
##
## Insert additional test code below but remember to change
## the print "1..x\n" in the BEGIN block to reflect the
## total number of tests that will be run. 

my $serror = 0;
my $ferror = 0;
my $format = '%-25s';

unless (eval "require SOAP::Lite; 1;") {
    print STDERR "SOAP::Lite not installed.\nThis means that Bio::Biblio module may not be usable. Skipping some tests.\n";
    $serror = 1;
}

use Bio::Root::IO;
my $testfile = Bio::Root::IO->catfile ('t','data','stress_test_medline.xml');
unless (-e $testfile) {
    print STDERR "Cannot find testing data '$testfile'. Skipping some tests.\n";
    $ferror = 1;
}

# check 'use ...'
eval { require Bio::Biblio };
print sprintf ($format, 'use Bio::Biblio'); skip ($error, %Bio::Biblio::);
print $@ if $@;

# check 'new...'
my $biblio;
print sprintf ($format, "new Bio::Biblio "); skip ($serror,
						   defined ($biblio = new Bio::Biblio (-location => 'http://localhost:4567')));

# check MEDLINE XML parser
eval { require Bio::Biblio::IO };
print sprintf ($format, "use Bio::Biblio::IO "); skip ($error, %Bio::Biblio::IO::);

my $io;
print sprintf ($format, "new Bio::Biblio::IO ");
skip ($ferror,
      defined (eval { $io = new Bio::Biblio::IO ('-format' => 'medlinexml',
						 '-file'   => $testfile,
						 '-result' => 'raw') }));
print "Reading and parsing XML file...\n";
print sprintf ($format, "    citation 1 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text1');
print sprintf ($format, "    citation 2 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text248');
print sprintf ($format, "    citation 3 "); skip ($ferror, eval { $io->next_bibref->{'medlineID'} }, 'Text495');

print "Getting citations using callback...\n";
my (@ids) = ('Text1', 'Text248', 'Text495');
my $callback_used = 'no';
unless ($ferror) {
    $io = new Bio::Biblio::IO ('-format'   => 'medlinexml',
			       '-file'     => $testfile,
#			       '-result'   => 'medline2ref',  # this is default
			       '-callback' => \&callback);
}
print sprintf ($format, "    calling callback "); skip ($ferror, $callback_used, 'yes');

sub callback {
    my $citation = shift;
    $callback_used = 'yes';
    print sprintf ($format, '    citation ' . (@ids+0) . ' '); skip ($ferror, $citation->{'_identifier'}, shift @ids);
}

print "Reading and parsing XML string...\n";
$io = new Bio::Biblio::IO ('-format'   => 'medlinexml',
			   '-data'     => <<XMLDATA,
<MedlineCitationSet>
<MedlineCitation>
<MedlineID>12345678</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
<MedlineCitation>
<MedlineID>abcdefgh</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
</MedlineCitationSet>
XMLDATA
			   '-result'   => 'medline2ref',
			   );
print sprintf ($format, "    citation 1 "); ok ($io->next_bibref->{'_identifier'}, '12345678');
print sprintf ($format, "    citation 2 "); ok ($io->next_bibref->{'_identifier'}, 'abcdefgh');

print "Reading and parsing XML string handle...\n";
use IO::String;
my $data = <<XMLDATA;
<MedlineCitationSet>
<MedlineCitation>
<MedlineID>87654321</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
<MedlineCitation>
<MedlineID>hgfedcba</MedlineID>
<Article><Journal></Journal></Article>
</MedlineCitation>
</MedlineCitationSet>
XMLDATA

$io = new Bio::Biblio::IO ('-format' => 'medlinexml',
			   '-fh'     => IO::String->new ($data),
			   );
print sprintf ($format, "    citation 1 "); ok ($io->next_bibref->{'_identifier'}, '87654321');
print sprintf ($format, "    citation 2 "); ok ($io->next_bibref->{'_identifier'}, 'hgfedcba');

if ($serror) {
    print "Remaining tests skipped.\n";
} else {

    print "Testing SOAP services...\n";

    # --- launch a testing SOAP server
    my ($pid, $port, $max_port);
    $port = 4444;
    $max_port = $port + 100;
    if ($pid = fork) {
	# parent here

	sleep 1;
	$biblio = new Bio::Biblio (-location => "tcp://localhost:$port", -namespace => 'soap_server');

	print sprintf ($format, '    get_count'); ok ($biblio->get_count, '43');
	print sprintf ($format, '    get_by_id'); ok ($biblio->get_by_id ('X'), 'X');
	print sprintf ($format, '    find (1)'); ok ($biblio->find ('a,b','c,d')->get_collection_id, 'a,b,c,d');
	print sprintf ($format, '    find (2)'); ok ($biblio->find (['x', 'y'], ['u', 'v'])->get_collection_id, 'x,y,u,v');

	print sprintf ($format, '    get_all_ids');
	ok ( eval { join (',', @{ $biblio->find ('AAA')->get_all_ids }) }, 'AAA'); print STDERR $@ if $@;
        
	print sprintf ($format, '    get_all');
	ok ( eval { join (',', @{ $biblio->find ('XXX')->get_all }) }, 'XXX'); print STDERR $@ if $@;
        
	print sprintf ($format, '    has_next');
	ok ( eval { $biblio->find (46)->has_next }, 1); print STDERR $@ if $@;
        
	print sprintf ($format, '    get_next');
	ok ( eval { $biblio->find ('BBB')->get_next }, 'BBB'); print STDERR $@ if $@;
        
	print sprintf ($format, '    get_more');
	ok ( eval { join (',', @{ $biblio->find ('CCC')->get_more (3) }) }, 'CCC,CCC,CCC'); print STDERR $@ if $@;
        
	print sprintf ($format, '    exists');
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
        sub getNext { shift; return shift; }
        sub getMore {
	    my ($self, $id, $how_many) = @_;
            my @result;
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
}

__END__
